library(SingleCellExperiment)
library(BPCells)
library(AnnotationDbi)
suppressPackageStartupMessages({
  library(dplyr)
  library(Rtsne)
  library(Matrix)
  library(uwot)
})

# building a pipeline roughly from here
# https://bnprks.github.io/BPCells/articles/pbmc3k.html#rna-normalization-pca-and-umap

run_bpcells <- function(
  sce, n_cluster, n_comp = 50, n_neig = 15, n_hvg = 1000,
  filter = c("manual", "auto"), time, clustering_info) {

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
  write(paste0("cells before: ", ncol(sce)), stderr())
  sce <- sce[, keep]
  write(paste0("cells after: ", ncol(sce)), stderr())
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Filter data. Time Elapsed:", time_elapsed))
  time$filter <- time_elapsed

  # NEED TO TIME THIS PART
  # sce is your existing SingleCellExperiment
  mat <- assay(sce, "counts")

  write_matrix_dir(mat = mat, dir = "bpcells_counts")
  bp_mat <- open_matrix_dir("bpcells_counts")

  # normalization ####
  start_time <- Sys.time()
  cell_sums <- Matrix::colSums(mat)
  mat_log <- mat %>%
    multiply_cols(1 / cell_sums) %>%
    `*`(10000) %>%
    log1p()
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Normalization. Time Elapsed:", time_elapsed))
  time$normalization <- time_elapsed

  # Identification of highly variable features (feature selection) ####
  start_time <- Sys.time()
  stats <- matrix_stats(mat_log, row_stats = "variance")
  gene_means <- stats$row_stats["mean", ]
  gene_vars  <- stats$row_stats["variance", ]

  # Vanilla choice: top HVGs by variance
  n_hvg <- min(n_hvg, nrow(mat_log))
  hvg_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_hvg)]
  hvg_idx <- sort(hvg_idx)

  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("Identification of highly variable features. Time Elapsed:", time_elapsed))
  time$hvg <- time_elapsed

  # PCA ####
  start_time <- Sys.time()
  mat_hvg <- mat_log[hvg_idx, ]
  mat_hvg <- write_matrix_memory(mat_hvg, compress = FALSE)
  svd <- BPCells::svds(mat_scaled, k = n_comp)

  # Cell embeddings: cells x PCs
  pca <- multiply_cols(svd$v, svd$d)
  pca <- as.matrix(pca)

  colnames(pca) <- paste0("PC_", seq_len(ncol(pca)))
  rownames(pca) <- colnames(mat)
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("PCA.Time Elapsed:", time_elapsed))
  time$pca <- time_elapsed

  # t-sne ####
  start_time <- Sys.time()
  set.seed(100000)
  tsne <- Rtsne::Rtsne(
    pca,
    dims = 2,
    perplexity = 30,
    check_duplicates = FALSE,
    pca = FALSE,
    verbose = TRUE
  )$Y
  rownames(tsne) <- rownames(pca)
  colnames(tsne) <- c("tSNE_1", "tSNE_2")
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("t-sne. Time Elapsed:", time_elapsed))
  time$t_sne <- time_elapsed

  # umap ####
  start_time <- Sys.time()
  set.seed(1000000)
  umap <- uwot::umap(
    pca,
    n_neighbors = 30,
    min_dist = 0.3,
    metric = "cosine",
    verbose = TRUE
  )
  rownames(umap) <- rownames(pca)
  colnames(umap) <- c("UMAP_1", "UMAP_2")
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste("umap. Time Elapsed:", time_elapsed))
  time$umap <- time_elapsed



# ----------------------------
# 7. Clustering
# ----------------------------
clust_louvain <- snn %>% cluster_graph_louvain()
clust_leiden  <- snn %>% cluster_graph_leiden()

  # louvain  ####
  start_time <- Sys.time()
  louvain_search <- binary_search(
    sce,
    do_clustering = function(snn, resolution) {
      clust_louvain <- snn %>% cluster_graph_louvain(resolution = resolution)
    },
    extract_nclust = function(result) length(unique(result)),
    n_clust_target = n_cluster
  )
  louvain_clustering <- louvain_search$result
  clustering_info$resolutions$louvain <- louvain_search$resolution
  clustering_info$num_runs$louvain <- louvain_search$num_runs
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  avg_time_elapsed <- time_elapsed / louvain_search$num_runs
  print(paste("Louvain clusterings. Total Search Time Elapsed:", time_elapsed))
  print(paste("Louvain clusterings. Average Time per Run:", avg_time_elapsed))
  print(paste("Louvain resolution:", clustering_info$resolutions$louvain))
  print(paste("Louvain runs:", clustering_info$num_runs$louvain))
  time$louvain <- avg_time_elapsed

  # leiden ####
  start_time <- Sys.time()
  knn <- knn_hnsw(pca, k = n_neig, ef = 200)
  snn <- knn %>% knn_to_snn_graph()
  leiden_search <- binary_search(
    sce,
    do_clustering = function(snn, resolution) {
      clust_louvain <- snn %>% cluster_graph_louvain(resolution = resolution)
    },
    extract_nclust = function(result) length(unique(result)),
    n_clust_target = n_cluster
  )
  leiden_clustering <- leiden_search$result
  clustering_info$resolutions$leiden <- leiden_search$resolution
  clustering_info$num_runs$leiden <- leiden_search$num_runs
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  avg_time_elapsed <- time_elapsed / leiden_search$num_runs
  print(paste("Leiden clusterings. Total Search Time Elapsed:", time_elapsed))
  print(paste("Leiden clusterings. Average Time per Run:", avg_time_elapsed))
  print(paste("Leiden resolution:", clustering_info$resolutions$leiden))
  print(paste("Leiden runs:", clustering_info$num_runs$leiden))
  time$leiden <- avg_time_elapsed

  return(list(
    pca = pca,
    hvgs = rownames(sce)[hvg_idx],
    cell_ids = colnames(sce),
    time = time,
    clustering_info = clustering_info,
    leiden = leiden_clustering,
    louvain = louvain_clustering
  ))
}
