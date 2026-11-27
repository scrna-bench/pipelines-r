#!/usr/bin/env Rscript
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "..", "src", "common.R"), chdir = TRUE)

# scrapper pipeline.

suppressPackageStartupMessages({
  library(DelayedArray)
  library(scrapper)
})

args <- pipeline_args()
sce <- read_input(args)
n_cluster <- n_cluster_target(args)
nthreads <- 1
assay(sce) <- DelayedArray(assay(sce))

qc_metrics <- with_phase("find_mit_gene", {
  computeRnaQcMetrics(assay(sce), subsets = list(mt = grepl("^MT-", rownames(sce))), num.threads = nthreads)
})

write(paste0("cells before: ", ncol(sce)), stderr())
keep <- with_phase("filter", {
  if (args$filter == "manual") {
    qc <- metadata(sce)$qc_thresholds
    mt_percent <- qc_metrics$subsets$mt * 100
    qc_metrics$detected >= qc[qc$metric == "nFeature", "min"] &
      qc_metrics$detected <= qc[qc$metric == "nFeature", "max"] &
      mt_percent < qc[qc$metric == "percent.mt", "max"] &
      qc_metrics$sum <= qc[qc$metric == "nCount", "max"]
  } else {
    filterRnaQcMetrics(suggestRnaQcThresholds(qc_metrics), qc_metrics)
  }
})
filtered <- sce[, keep, drop = FALSE]
write(paste0("cells after: ", ncol(filtered)), stderr())

assay(filtered, "normalized") <- with_phase("normalization", {
  normalizeCounts(assay(filtered), centerSizeFactors(qc_metrics$sum[keep]))
})

hvgs <- with_phase("hvg", {
  gene_var <- modelGeneVariances(assay(filtered, "normalized"), num.threads = nthreads)
  chooseHighlyVariableGenes(gene_var$statistics$residuals, top = args$n_hvg)
})

pca <- with_phase("pca", {
  runPca(assay(filtered, "normalized")[hvgs, ], num.threads = nthreads, number = args$n_comp)
})
reducedDim(filtered, "PCA") <- t(pca$components)

reducedDim(filtered, "TSNE") <- with_phase("t_sne", runTsne(pca$components, num.threads = nthreads))
reducedDim(filtered, "UMAP") <- with_phase("umap", {
  set.seed(1000000)
  runUmap(pca$components, num.threads = nthreads)
})

louvain <- search_phase("louvain", {
  graph <- buildSnnGraph(pca$components, num.neighbors = args$n_neig, num.threads = nthreads)
  binary_search(
    graph,
    do_clustering = function(spe, resolution) {
      clusterGraph(spe, method = "multilevel", multilevel.resolution = resolution)
    },
    extract_nclust = function(result) length(unique(result$membership)),
    n_clust_target = n_cluster
  )
})
leiden <- search_phase("leiden", {
  graph <- buildSnnGraph(pca$components, num.neighbors = args$n_neig, num.threads = nthreads)
  binary_search(
    graph,
    do_clustering = function(spe, resolution) {
      clusterGraph(spe, method = "leiden", leiden.resolution = resolution)
    },
    extract_nclust = function(result) length(unique(result$membership)),
    n_clust_target = n_cluster
  )
})

write_outputs(
  args,
  cell_ids = colnames(filtered),
  pca = reducedDim(filtered, "PCA"),
  hvgs = rownames(filtered)[hvgs],
  louvain = louvain$result$membership,
  leiden = leiden$result$membership,
  info = clustering_info(louvain, leiden)
)
