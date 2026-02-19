library(SingleCellExperiment)
library(scuttle)
library(AnnotationDbi)
library(scran)
library(scater)
library(bluster)
library(EnsDb.Hsapiens.v75)

run_osca <- function(
  sce, n_cluster, n_comp = 50, n_neig = 15, n_hvg = 1000,
  filter = c("manual", "auto"), time, resolutions
) {
  filter <- match.arg(filter)
  #### 1. find mithocondial genes  ####
  start_time <- Sys.time()
  chr.loc <- mapIds(EnsDb.Hsapiens.v75,
    keys = rownames(sce),
    keytype = "SYMBOL", column = "SEQNAME"
  )
  is.mito <- which(chr.loc == "MT")
  df <- perCellQCMetrics(sce, subsets = list(Mito = is.mito))
  # include them in the object
  colData(sce) <- cbind(colData(sce), df)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Find mithocondial genes. Time Elapsed:", time_elapsed))
  time$find_mit_gene <- time_elapsed

  # 2. filter data ####
  start_time <- Sys.time()
  if (filter == "manual") {
    qc <- metadata(sce)$qc_thresholds
    keep <- df$detected > qc[qc$metric == "nFeature", "min"] &
      df$detected < qc[qc$metric == "nFeature", "max"] &
      df$subsets_Mito_percent < qc[qc$metric == "percent.mt", "max"] &
      df$sum < qc[qc$metric == "nCount", "max"]
  } else {
    reasons <- perCellQCFilters(df, sub.fields = "subsets_Mito_percent")
    keep <- !reasons$discard
  }
  sce <- sce[, keep]
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Filter data. Time Elapsed:", time_elapsed))
  time$filter <- time_elapsed

  # normalization ####
  start_time <- Sys.time()
  sce <- logNormCounts(sce)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Normalization. Time Elapsed:", time_elapsed))
  time$normalization <- time_elapsed

  # Identification of highly variable features (feature selection) ####
  start_time <- Sys.time()
  dec.sce <- modelGeneVar(sce)
  hvg.sce.var <- getTopHVGs(dec.sce, n = n_hvg)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Identification of highly variable features. Time Elapsed:", time_elapsed))
  time$hvg <- time_elapsed

  # PCA ####
  start_time <- Sys.time()
  sce <- runPCA(sce, subset_row = hvg.sce.var, ncomponents = n_comp)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("PCA.Time Elapsed:", time_elapsed))
  time$pca <- time_elapsed

  # t-sne ####
  start_time <- Sys.time()
  set.seed(100000)
  sce <- runTSNE(sce, dimred = "PCA")
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("t-sne. Time Elapsed:", time_elapsed))
  time$t_sne <- time_elapsed

  # umap ####
  start_time <- Sys.time()
  set.seed(1000000)
  sce <- runUMAP(sce, dimred = "PCA")
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("umap. Time Elapsed:", time_elapsed))
  time$umap <- time_elapsed

  # louvain  ####
  start_time <- Sys.time()
  louvain_search <- binary_search(
    sce,
    do_clustering = function(spe, resolution) {
      clusterCells(spe,
        use.dimred = "PCA",
        BLUSPARAM = NNGraphParam(
          k = n_neig,
          cluster.fun = "louvain", cluster.args = list(resolution = resolution)
        )
      )
    },
    extract_nclust = function(result) length(unique(result)),
    n_clust_target = n_cluster
  )
  louvain_clustering <- louvain_search$result
  resolutions$louvain <- louvain_search$resolution
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Louvain clusterings. Time Elapsed:", time_elapsed))
  print(paste("Louvain resolution:", resolutions$louvain))
  time$louvain <- time_elapsed

  # leiden ####
  start_time <- Sys.time()
  leiden_search <- binary_search(
    sce,
    do_clustering = function(spe, resolution) {
      clusterCells(spe,
        use.dimred = "PCA",
        BLUSPARAM = NNGraphParam(
          k = n_neig,
          cluster.fun = "leiden", cluster.args = list(resolution = resolution)
        )
      )
    },
    extract_nclust = function(result) length(unique(result)),
    n_clust_target = n_cluster
  )
  leiden_clustering <- leiden_search$result
  resolutions$leiden <- leiden_search$resolution
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Leiden clusterings. Time Elapsed:", time_elapsed))
  print(paste("Leiden resolution:", resolutions$leiden))
  time$leiden <- time_elapsed

  return(list(
    pca = reducedDim(sce, "PCA"),
    hvgs = hvg.sce.var,
    cell_ids = colnames(sce),
    time = time,
    resolutions = resolutions,
    leiden = leiden_clustering,
    louvain = louvain_clustering
  ))
}
