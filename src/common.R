# Shared by every pipeline entrypoint: arguments, input loading, the target
# number of clusters, and writing the stage outputs. Source with chdir = TRUE.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

source("cli.R")
source("timing.R")
source("search.R")

# Parse the CLI (inputs, parameters, outputs), start the event log, and return
# a list with typed parameters and the output paths.
pipeline_args <- function(filters = c("manual", "auto")) {
  a <- parse_cli()
  require_args(a, c("data_h5ad", "data_clusters_truth_num", "d_cluster"))
  args <- list(
    data_h5ad = a$data_h5ad,
    data_clusters_truth_num = a$data_clusters_truth_num,
    d_cluster = as.integer(a$d_cluster),
    n_comp = as.integer(a$n_comp %||% 50),
    n_neig = as.integer(a$n_neig %||% 15),
    n_hvg = as.integer(a$n_hvg %||% 1000),
    filter = match.arg(a$filter %||% filters[1], filters),
    out = list(
      clusters = output_path(a, "methods_clusters", ".clusters.tsv"),
      pca = output_path(a, "methods_pca", ".pca.tsv"),
      hvgs = output_path(a, "methods_hvgs", ".hvgs.tsv"),
      clustering_info = output_path(a, "methods_clustering_info", ".clustering_info.json"),
      events = a$outputs$methods_events %||% file.path(a$output_dir, "obkit-events.jsonl")
    )
  )
  timing_init(args$out$events)
  args
}

# The input h5ad as a SingleCellExperiment with a `counts` assay.
read_input <- function(args) {
  with_phase("load", {
    sce <- anndataR::read_h5ad(args$data_h5ad, as = "SingleCellExperiment")
    names(assays(sce))[names(assays(sce)) == "X"] <- "counts"
    sce
  })
}

n_cluster_target <- function(args) {
  as.integer(readLines(args$data_clusters_truth_num, n = 1)) + args$d_cluster
}

# Run a binary_search() inside a timed phase; the end event carries the chosen
# resolution and the number of clustering runs.
search_phase <- function(name, search) {
  with_phase(name, search, attrs = function(s) {
    list(resolution = s$resolution, num_runs = s$num_runs)
  })
}

clustering_info <- function(louvain, leiden) {
  list(
    resolutions = list(louvain = louvain$resolution, leiden = leiden$resolution),
    num_runs = list(louvain = louvain$num_runs, leiden = leiden$num_runs)
  )
}

# pca: cells x PCs matrix with cell ids as rownames.
write_outputs <- function(args, cell_ids, pca, hvgs, louvain, leiden, info) {
  with_phase("write", {
    write.table(
      data.frame(cell_id = cell_ids, louvain = louvain, leiden = leiden),
      args$out$clusters,
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    pca <- as.matrix(pca)
    colnames(pca) <- paste0("PC", seq_len(ncol(pca)))
    write.table(
      data.frame(cell_id = rownames(pca), pca, check.names = FALSE),
      args$out$pca,
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    writeLines(as.character(hvgs), args$out$hvgs)
    jsonlite::write_json(info, args$out$clustering_info, auto_unbox = TRUE, pretty = TRUE)
  })
}
