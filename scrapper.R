library(SingleCellExperiment)
library(DelayedArray)
library(scrapper)

run_scrapper <- function(
  sce, n_cluster, n_comp = 50, n_neig = 15, n_hvg = 1000,
  filter = c("manual", "auto"), time, resolutions
) {
  filter <- match.arg(filter)
  nthreads <- 1
  assay(sce) <- DelayedArray(assay(sce))

  #### 1. find mithocondial genes  ####
  start_time <- Sys.time()
  is.mito <- grepl("^MT-", rownames(sce))
  rna.qc.metrics <- computeRnaQcMetrics(assay(sce),
    subsets = list(mt = is.mito),
    num.threads = nthreads
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Find mithocondrial genes. Time Elapsed:", time_elapsed))
  time$find_mit_gene <- time_elapsed

  # filter data ####
  start_time <- Sys.time()
  if (filter == "manual") {
    qc <- metadata(sce)$qc_thresholds
    mt_percent <- rna.qc.metrics$subsets$mt * 100
    keep <- rna.qc.metrics$detected > qc[qc$metric == "nFeature", "min"] &
      rna.qc.metrics$detected < qc[qc$metric == "nFeature", "max"] &
      mt_percent < qc[qc$metric == "percent.mt", "max"] &
      rna.qc.metrics$sum < qc[qc$metric == "nCount", "max"]
  } else {
    rna.qc.thresholds <- suggestRnaQcThresholds(rna.qc.metrics)
    keep <- filterRnaQcMetrics(rna.qc.thresholds, rna.qc.metrics)
  }
  filtered <- sce[, keep, drop = FALSE]
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Filter data. Time Elapsed:", time_elapsed))
  time$filter <- time_elapsed

  # normalization ####
  start_time <- Sys.time()
  size.factors <- centerSizeFactors(rna.qc.metrics$sum[keep])
  normalized <- normalizeCounts(assay(filtered), size.factors)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Normalization. Time Elapsed:", time_elapsed))
  time$normalization <- time_elapsed
  assay(filtered, "normalized") <- normalized

  # Identification of highly variable features (feature selection) ####
  start_time <- Sys.time()
  gene.var <- modelGeneVariances(
    assay(filtered, "normalized"),
    num.threads = nthreads
  )
  hvg.sce.var <- chooseHighlyVariableGenes(
    gene.var$statistics$residuals,
    top = n_hvg
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Identification of highly variable features. Time Elapsed:", time_elapsed))
  time$hvg <- time_elapsed

  # PCA ####
  start_time <- Sys.time()
  pca <- runPca(
    (assay(filtered, "normalized")[hvg.sce.var, ]),
    num.threads = nthreads, number = n_comp
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("PCA. Time Elapsed:", time_elapsed))
  time$pca <- time_elapsed
  reducedDim(filtered, "PCA") <- t(pca$components)

  # t-sne ####
  start_time <- Sys.time()
  tsne.out <- runTsne(pca$components, num.threads = nthreads)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("t-sne. Time Elapsed:", time_elapsed))
  time$t_sne <- time_elapsed
  reducedDim(filtered, "TSNE") <- tsne.out

  # umap ####
  start_time <- Sys.time()
  set.seed(1000000)
  umap.out <- runUmap(pca$components, num.threads = nthreads)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("UMAP. Time Elapsed:", time_elapsed))
  time$umap <- time_elapsed
  reducedDim(filtered, "UMAP") <- umap.out

  # louvain  ####
  start_time <- Sys.time()
  snn.graph.louvain <- buildSnnGraph(
    pca$components,
    num.neighbors = n_neig, num.threads = nthreads
  )
  louvain_search <- binary_search(
    snn.graph.louvain,
    do_clustering = function(spe, resolution) {
      clusterGraph(spe, method = c("multilevel"), multilevel.resolution = resolution)
    },
    extract_nclust = function(result) length(unique(result$membership)),
    n_clust_target = n_cluster
  )
  louvain_clustering <- louvain_search$result$membership
  resolutions$louvain <- louvain_search$resolution
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Louvain clusterings. Time Elapsed:", time_elapsed))
  print(paste("Louvain resolution:", resolutions$louvain))
  time$louvain <- time_elapsed

  # leiden ####
  start_time <- Sys.time()
  snn.graph.leiden <- buildSnnGraph(
    pca$components,
    num.neighbors = n_neig, num.threads = nthreads
  )
  leiden_search <- binary_search(
    snn.graph.leiden,
    do_clustering = function(spe, resolution) {
      clusterGraph(spe, method = c("leiden"), leiden.resolution = resolution)
    },
    extract_nclust = function(result) length(unique(result$membership)),
    n_clust_target = n_cluster
  )
  leiden_clustering <- leiden_search$result$membership
  resolutions$leiden <- leiden_search$resolution
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Leiden clusterings. Time Elapsed:", time_elapsed))
  print(paste("Leiden resolution:", resolutions$leiden))
  time$leiden <- time_elapsed

  return(list(
    pca = reducedDim(filtered, "PCA"),
    hvgs = rownames(filtered)[hvg.sce.var],
    cell_ids = colnames(filtered),
    time = time,
    resolutions = resolutions,
    leiden = leiden_clustering,
    louvain = louvain_clustering
  ))
}
