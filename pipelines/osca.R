#!/usr/bin/env Rscript
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "..", "src", "common.R"), chdir = TRUE)

# OSCA (scran/scater/bluster) pipeline.

suppressPackageStartupMessages({
  library(scuttle)
  library(AnnotationDbi)
  library(scran)
  library(scater)
  library(bluster)
  library(EnsDb.Hsapiens.v75)
})

args <- pipeline_args()
sce <- read_input(args)
n_cluster <- n_cluster_target(args)

df <- with_phase("find_mit_gene", {
  chr_loc <- mapIds(EnsDb.Hsapiens.v75, keys = rownames(sce), keytype = "SYMBOL", column = "SEQNAME")
  df <- perCellQCMetrics(sce, subsets = list(Mito = which(chr_loc == "MT")))
  colData(sce) <- cbind(colData(sce), df)
  df
})

write(paste0("cells before: ", ncol(sce)), stderr())
sce <- with_phase("filter", {
  if (args$filter == "manual") {
    qc <- metadata(sce)$qc_thresholds
    keep <- df$detected > qc[qc$metric == "nFeature", "min"] &
      df$detected < qc[qc$metric == "nFeature", "max"] &
      df$subsets_Mito_percent < qc[qc$metric == "percent.mt", "max"] &
      df$sum < qc[qc$metric == "nCount", "max"]
  } else {
    keep <- !perCellQCFilters(df, sub.fields = "subsets_Mito_percent")$discard
  }
  sce[, keep]
})
write(paste0("cells after: ", ncol(sce)), stderr())

sce <- with_phase("normalization", logNormCounts(sce))
hvgs <- with_phase("hvg", getTopHVGs(modelGeneVar(sce), n = args$n_hvg))
sce <- with_phase("pca", runPCA(sce, subset_row = hvgs, ncomponents = args$n_comp))
sce <- with_phase("t_sne", {
  set.seed(100000)
  runTSNE(sce, dimred = "PCA")
})
sce <- with_phase("umap", {
  set.seed(1000000)
  runUMAP(sce, dimred = "PCA")
})

cluster_fn <- function(algorithm) {
  function(spe, resolution) {
    clusterCells(spe,
      use.dimred = "PCA",
      BLUSPARAM = NNGraphParam(
        k = args$n_neig,
        cluster.fun = algorithm, cluster.args = list(resolution = resolution)
      )
    )
  }
}
louvain <- search_phase("louvain", binary_search(
  sce, cluster_fn("louvain"), function(result) length(unique(result)), n_cluster
))
leiden <- search_phase("leiden", binary_search(
  sce, cluster_fn("leiden"), function(result) length(unique(result)), n_cluster
))

write_outputs(
  args,
  cell_ids = colnames(sce),
  pca = reducedDim(sce, "PCA"),
  hvgs = hvgs,
  louvain = louvain$result,
  leiden = leiden$result,
  info = clustering_info(louvain, leiden)
)
