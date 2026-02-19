#!/usr/bin/env Rscript

library(argparse)
library(jsonlite)
library(zellkonverter)

parser <- ArgumentParser(description = "Benchmarking entrypoint")

parser$add_argument(
  "--output_dir", "-o",
  dest = "output_dir", type = "character",
  help = "output directory where files will be saved",
  default = getwd(), required = TRUE
)
parser$add_argument(
  "--name", "-n",
  dest = "name", type = "character",
  help = "name of the module",
  required = TRUE
)
parser$add_argument(
  "--data.h5ad",
  dest = "data_path", type = "character",
  help = "input data h5ad path",
  required = TRUE
)
parser$add_argument(
  "--method_name",
  dest = "method_name", type = "character",
  help = "name of the method",
  choices = c("osca", "scrapper", "seurat"), required = TRUE
)
parser$add_argument(
  "--n_cluster",
  dest = "n_cluster", type = "integer",
  help = "target number of clusters",
  required = TRUE
)
parser$add_argument(
  "--n_comp",
  dest = "n_comp", type = "integer",
  help = "number of PCA components to use for KNN graph construction",
  default = 50, required = FALSE
)
parser$add_argument(
  "--n_neig",
  dest = "n_neig", type = "integer",
  help = "number of neighbors to use for KNN graph construction",
  default = 15, required = FALSE
)
parser$add_argument(
  "--n_hvg",
  dest = "n_hvg", type = "integer",
  help = "number of highly variable genes to use",
  default = 1000, required = FALSE
)
parser$add_argument(
  "--filter",
  dest = "filter", type = "character",
  help = "filtering (manual uses suggested cutoffs, auto uses package QC)",
  choices = c("manual", "auto"),
  default = "manual", required = FALSE
)

args <- parser$parse_args()

cargs <- commandArgs(trailingOnly = FALSE)
m <- grep("--file=", cargs)
run_dir <- dirname(gsub("--file=", "", cargs[[m]]))

timings_path <- file.path(args$output_dir, paste0(args$name, ".timings.json"))
resolutions_path <- file.path(args$output_dir, paste0(args$name, ".resolutions.json"))
clusters_path <- file.path(args$output_dir, paste0(args$name, ".clusters.tsv"))
pca_path <- file.path(args$output_dir, paste0(args$name, ".pca.tsv"))
hvgs_path <- file.path(args$output_dir, paste0(args$name, ".hvgs.tsv"))

search_res_path <- file.path(run_dir, "search_res.r")
source(search_res_path)

sce <- readH5AD(args$data_path, reader = "python")

# time object to store time involved (in seconds) in each step
time <- list(
  find_mit_gene = NA_real_, filter = NA_real_, normalization = NA_real_,
  hvg = NA_real_, scaling = NA_real_, pca = NA_real_,
  t_sne = NA_real_, umap = NA_real_,
  louvain = NA_real_, leiden = NA_real_
)

resolutions <- list(louvain = NA_real_, leiden = NA_real_)

# source and run appropriate method
if (args$method_name == "seurat") {
  seurat_r_path <- file.path(run_dir, "seurat.R")
  source(seurat_r_path)
  output_data <- run_seurat(
    sce,
    args$n_cluster, args$n_comp, args$n_neig, args$n_hvg, args$filter, time, resolutions
  )
} else if (args$method_name == "osca") {
  osca_r_path <- file.path(run_dir, "OSCA.R")
  source(osca_r_path)
  output_data <- run_osca(
    sce,
    args$n_cluster, args$n_comp, args$n_neig, args$n_hvg, args$filter, time, resolutions
  )
} else if (args$method_name == "scrapper") {
  scrapper_r_path <- file.path(run_dir, "scrapper.R")
  source(scrapper_r_path)
  output_data <- run_scrapper(
    sce,
    args$n_cluster, args$n_comp, args$n_neig, args$n_hvg, args$filter, time, resolutions
  )
}

# write outputs to files
output_data$time <- lapply(output_data$time, function(x) {
  as.numeric(x, units = "secs")
})
write_json(
  output_data$time, timings_path,
  auto_unbox = TRUE, pretty = TRUE
)
write_json(
  output_data$resolutions, resolutions_path,
  auto_unbox = TRUE, pretty = TRUE
)
write.table(
  data.frame(
    cell_id = output_data$cell_ids,
    louvain = output_data$louvain, leiden = output_data$leiden
  ),
  clusters_path,
  sep = "\t", quote = F, row.names = F
)
output_data$pca <- data.frame(
  cell_id = rownames(output_data$pca),
  output_data$pca
)
colnames(output_data$pca) <- c(
  "cell_id",
  paste0("PC", seq_len(ncol(output_data$pca) - 1))
)
write.table(
  output_data$pca, pca_path,
  sep = "\t", quote = F, row.names = F
)
cat(str(output_data))
cat(str(output_data$hvgs))
writeLines(output_data$hvgs, file(hvgs_path))
