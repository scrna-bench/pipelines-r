#!/usr/bin/env Rscript
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "..", "src", "common.R"), chdir = TRUE)

# bixverse pipeline (Rust-backed). Ported by Mark Robinson in
# omni-training/pipelines-r (bixverse.R).

suppressPackageStartupMessages({
  library(bixverse)
})

args <- pipeline_args()
sce <- read_input(args)
n_cluster <- n_cluster_target(args)

# caps bixverse's Rust thread pool if it has not started yet in this session
Sys.setenv(RAYON_NUM_THREADS = 1)

qc_thresholds <- metadata(sce)$qc_thresholds
colData(sce)$cell_id <- colnames(sce)
rowData(sce)$gene_id <- rownames(sce)
dir_data <- tempfile("bixverse_")
dir.create(dir_data, recursive = TRUE)

# bixverse normalises while loading the counts into its on-disk store, so this
# phase stands in for normalisation
sc_object <- with_phase("normalization", load_sce(
  object = SingleCells(dir_data = dir_data),
  sce = sce,
  assay_name = "counts",
  # minimal thresholds so the loader keeps (almost) everything; the real
  # filtering happens below
  sc_qc_param = params_sc_min_quality(min_unique_genes = 1L, min_lib_size = 1L, min_cells = 1L),
  .verbose = FALSE
))

sc_object <- with_phase("find_mit_gene", {
  var <- get_sc_var(sc_object)
  gene_set_proportions_sc(
    sc_object,
    list(MT = var[grepl("^MT-", gene_id), gene_id]),
    streaming = FALSE,
    .verbose = FALSE
  )
})

sc_object <- with_phase("filter", {
  qc_df <- sc_object[[c("cell_id", "lib_size", "nnz", "MT")]]
  write(paste0("cells before: ", nrow(qc_df)), stderr())
  if (args$filter == "manual") {
    qc <- qc_thresholds
    keep <- qc_df$nnz >= qc[qc$metric == "nFeature", "min"] &
      qc_df$nnz <= qc[qc$metric == "nFeature", "max"] &
      qc_df$MT * 100 < qc[qc$metric == "percent.mt", "max"] &
      qc_df$lib_size <= qc[qc$metric == "nCount", "max"]
  } else {
    # MAD-based outliers, as in the bixverse PBMC vignette
    cell_qc <- run_cell_qc(
      metrics = list(
        log10_lib_size = log10(qc_df$lib_size),
        log10_nnz = log10(qc_df$nnz),
        MT = qc_df$MT
      ),
      cells_to_keep = get_cells_to_keep(sc_object),
      directions = c(log10_lib_size = "twosided", log10_nnz = "twosided", MT = "above"),
      threshold = 3
    )
    keep <- !cell_qc$combined
  }
  sc_object <- set_cells_to_keep(sc_object, qc_df[keep, cell_id])
  write(paste0("cells after: ", length(get_cells_to_keep(sc_object))), stderr())
  sc_object
})

sc_object <- with_phase("hvg", find_hvg_sc(object = sc_object, hvg_no = args$n_hvg, .verbose = FALSE))
sc_object <- with_phase("pca", calculate_pca_sc(object = sc_object, no_pcs = args$n_comp, sparse_svd = TRUE))
sc_object <- with_phase("t_sne", tsne_sc(sc_object, .verbose = FALSE))
# use_knn = FALSE: neighbours come from the PCA, as in runUMAP()
sc_object <- with_phase("umap", umap_sc(sc_object, use_knn = FALSE, seed = 42L, .verbose = FALSE))

cluster_fn <- function(algorithm) {
  function(spe, resolution) {
    spe <- find_clusters_sc(spe, cluster_algorithm = algorithm, res = resolution, name = algorithm)
    get_sc_obs(spe, filtered = TRUE)[[algorithm]]
  }
}
find_neighbours <- function(sc_object) {
  find_neighbours_sc(
    object = sc_object,
    neighbours_params = params_sc_neighbours(knn = list(k = args$n_neig)),
    .verbose = FALSE
  )
}
louvain <- search_phase("louvain", {
  sc_object <- find_neighbours(sc_object)
  binary_search(sc_object, cluster_fn("louvain"), function(result) length(unique(result)), n_cluster)
})
leiden <- search_phase("leiden", {
  sc_object <- find_neighbours(sc_object)
  binary_search(sc_object, cluster_fn("leiden"), function(result) length(unique(result)), n_cluster)
})

cell_ids <- get_sc_obs(sc_object, filtered = TRUE)[["cell_id"]]
pca <- get_pca_factors(sc_object)
rownames(pca) <- cell_ids

write_outputs(
  args,
  cell_ids = cell_ids,
  pca = pca,
  hvgs = unname(get_gene_names_from_idx(sc_object, get_hvg(sc_object))),
  louvain = louvain$result,
  leiden = leiden$result,
  info = clustering_info(louvain, leiden)
)
unlink(dir_data, recursive = TRUE, force = TRUE)
