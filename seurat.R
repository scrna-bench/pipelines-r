library(Seurat)
library(SingleCellExperiment)

run_seurat <- function(
  sce, n_cluster, n_comp = 50, n_neig = 15, n_hvg = 1000,
  filter = c("manual", "auto"), time, resolutions
) {
  filter <- match.arg(filter)
  # data ####
  data <- as.Seurat(sce, counts = "counts", data = NULL, assay = NULL)
  DefaultAssay(data) <- "originalexp"

  # find mitocondrial genes ####
  start_time <- Sys.time()
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Find mitocondrial genes. Time Elapsed:", time_elapsed))
  time$find_mit_gene <- time_elapsed

  # filter data ####
  write(paste0("before: ", dim(data)), stderr())
  start_time <- Sys.time()
  if (filter == "manual") {
    qc <- metadata(sce)$qc_thresholds
    data <- subset(
      data,
      subset = nFeature_originalexp > qc[qc$metric == "nFeature", "min"] &
        nFeature_originalexp < qc[qc$metric == "nFeature", "max"] &
        percent.mt < qc[qc$metric == "percent.mt", "max"] &
        nCount_originalexp < qc[qc$metric == "nCount", "max"]
    )
  } else {
    # seurat auto pipeline uses scuttle for filtering
    is.mito <- grepl("^MT-", rownames(sce))
    df <- scuttle::perCellQCMetrics(sce, subsets = list(Mito = is.mito))
    reasons <- scuttle::perCellQCFilters(
      df,
      sub.fields = "subsets_Mito_percent"
    )
    keep_cells <- colnames(sce)[!reasons$discard]
    data <- subset(data, cells = keep_cells)
  }
  end_time <- Sys.time()
  write(paste0("after: ", dim(data)), stderr())
  time_elapsed <- end_time - start_time
  print(paste("Filter data. Time Elapsed:", time_elapsed))
  time$filter <- time_elapsed

  # normalization ####
  start_time <- Sys.time()
  data <- NormalizeData(data, normalization.method = "LogNormalize")
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Normalization. Time Elapsed:", time_elapsed))
  time$normalization <- time_elapsed

  # Identification of highly variable features (feature selection) ####
  start_time <- Sys.time()
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = n_hvg)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Identification of highly variable features. Time Elapsed:", time_elapsed))
  time$hvg <- time_elapsed

  # Scaling the data ####
  start_time <- Sys.time()
  data <- ScaleData(data)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Scaling the data. Time Elapsed:", time_elapsed))
  time$scaling <- time_elapsed

  # PCA ####
  start_time <- Sys.time()
  data <- RunPCA(
    data,
    features = VariableFeatures(object = data), npcs = n_comp, verbose = FALSE
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("PCA. Time Elapsed:", time_elapsed))
  time$pca <- time_elapsed

  # t-sne ####
  start_time <- Sys.time()
  data <- RunTSNE(data,
    reduction = "pca", perplexity = 18
  )
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("t-sne. Time Elapsed:", time_elapsed))
  time$t_sne <- time_elapsed

  # UMAP ####
  start_time <- Sys.time()
  data <- RunUMAP(data, dims = 1:n_comp)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("UMAP. Time Elapsed:", time_elapsed))
  time$umap <- time_elapsed

  # louvain ####
  start_time <- Sys.time()
  data <- FindNeighbors(
    data,
    dims = 1:n_comp, k.param = n_neig, verbose = T
  )
  louvain_search <- binary_search(
    data,
    do_clustering = function(spe, resolution) {
      FindClusters(spe, algorithm = 1, cluster.name = "louvain", resolution = resolution)
    },
    extract_nclust = function(result) length(unique(result$louvain)),
    n_clust_target = n_cluster
  )
  data <- louvain_search$result
  resolutions$louvain <- louvain_search$resolution
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Louvain Clustering. Time Elapsed:", time_elapsed))
  print(paste("Louvain resolution:", resolutions$louvain))
  time$louvain <- time_elapsed

  # leiden ####
  start_time <- Sys.time()
  leiden_search <- binary_search(
    data,
    do_clustering = function(spe, resolution) {
      FindClusters(spe, algorithm = 4, cluster.name = "leiden", resolution = resolution)
    },
    extract_nclust = function(result) length(unique(result$leiden)),
    n_clust_target = n_cluster
  )
  data <- leiden_search$result
  resolutions$leiden <- leiden_search$resolution
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Leiden Clustering. Time Elapsed:", time_elapsed))
  print(paste("Leiden resolution:", resolutions$leiden))
  time$leiden <- time_elapsed

  return(list(
    pca = Embeddings(data, reduction = "pca"),
    hvgs = VariableFeatures(data),
    time = time,
    resolutions = resolutions,
    cell_ids = colnames(data),
    leiden = data$leiden,
    louvain = data$louvain
  ))
}
