#!/usr/bin/env Rscript
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "..", "src", "common.R"), chdir = TRUE)

# Seurat pipeline.

suppressPackageStartupMessages({
  library(Seurat)
})

args <- pipeline_args()
sce <- read_input(args)
n_cluster <- n_cluster_target(args)

data <- as.Seurat(sce, counts = "counts", data = NULL, assay = NULL)
DefaultAssay(data) <- "originalexp"

with_phase("find_mit_gene", {
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
})

write(paste0("cells before: ", ncol(data)), stderr())
data <- with_phase("filter", {
  if (args$filter == "manual") {
    qc <- metadata(sce)$qc_thresholds
    subset(
      data,
      subset = nFeature_originalexp > qc[qc$metric == "nFeature", "min"] &
        nFeature_originalexp < qc[qc$metric == "nFeature", "max"] &
        percent.mt < qc[qc$metric == "percent.mt", "max"] &
        nCount_originalexp < qc[qc$metric == "nCount", "max"]
    )
  } else {
    # the auto pipeline uses scuttle's outlier filters
    is_mito <- grepl("^MT-", rownames(sce))
    df <- scuttle::perCellQCMetrics(sce, subsets = list(Mito = is_mito))
    reasons <- scuttle::perCellQCFilters(df, sub.fields = "subsets_Mito_percent")
    subset(data, cells = colnames(sce)[!reasons$discard])
  }
})
write(paste0("cells after: ", ncol(data)), stderr())

data <- with_phase("normalization", NormalizeData(data, normalization.method = "LogNormalize"))
data <- with_phase("hvg", FindVariableFeatures(data, selection.method = "vst", nfeatures = args$n_hvg))
data <- with_phase("scaling", ScaleData(data))
data <- with_phase("pca", RunPCA(data, features = VariableFeatures(object = data), npcs = args$n_comp, verbose = FALSE))
data <- with_phase("t_sne", RunTSNE(data, reduction = "pca", perplexity = 18))
data <- with_phase("umap", RunUMAP(data, dims = 1:args$n_comp))

louvain <- search_phase("louvain", {
  data <- FindNeighbors(data, dims = 1:args$n_comp, k.param = args$n_neig, verbose = TRUE)
  binary_search(
    data,
    do_clustering = function(spe, resolution) {
      FindClusters(spe, algorithm = 1, cluster.name = "louvain", resolution = resolution)
    },
    extract_nclust = function(result) length(unique(result$louvain)),
    n_clust_target = n_cluster
  )
})
data <- louvain$result

leiden <- search_phase("leiden", binary_search(
  data,
  do_clustering = function(spe, resolution) {
    FindClusters(spe, algorithm = 4, cluster.name = "leiden", resolution = resolution)
  },
  extract_nclust = function(result) length(unique(result$leiden)),
  n_clust_target = n_cluster
))
data <- leiden$result

write_outputs(
  args,
  cell_ids = colnames(data),
  pca = Embeddings(data, reduction = "pca"),
  hvgs = VariableFeatures(data),
  louvain = data$louvain,
  leiden = data$leiden,
  info = clustering_info(louvain, leiden)
)
