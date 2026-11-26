# Command-line parsing for omnibenchmark (api 0.8) entrypoints, base R only.
#
# ob calls an entrypoint as
#   <script> --output_dir D --name N --output <id>=<path> ... --<input> <path> --<param> <value>
# parse_cli() returns a named list of every `--key value` / `--key=value` pair,
# plus `$outputs`, a named list output id -> path built from the repeated
# `--output id=PATH` flags. Values stay character: callers convert types.

parse_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
  args <- list()
  outputs <- list()
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--")) stop(sprintf("unexpected argument '%s'", a), call. = FALSE)
    if (grepl("=", a, fixed = TRUE)) {
      key <- sub("^--([^=]+)=.*$", "\\1", a)
      val <- sub("^--[^=]+=", "", a)
      i <- i + 1L
    } else {
      key <- sub("^--", "", a)
      val <- if (i < length(argv)) argv[[i + 1L]] else NA_character_
      i <- i + 2L
    }
    if (key == "output") {
      if (is.na(val) || !grepl("^[^=]+=.+", val)) {
        stop(sprintf("--output wants ID=PATH, got '%s'", val), call. = FALSE)
      }
      eq <- regexpr("=", val, fixed = TRUE)
      outputs[[substr(val, 1L, eq - 1L)]] <- substring(val, eq + 1L)
    } else {
      args[[key]] <- val
    }
  }
  args$outputs <- outputs
  args
}

# Stop unless every flag in `keys` was given.
require_args <- function(args, keys) {
  missing <- setdiff(keys, names(args))
  if (length(missing)) {
    stop("missing required argument(s): ", paste0("--", missing, collapse = ", "), call. = FALSE)
  }
  invisible(args)
}

# Path for output `id`: `--output id=PATH` when given, else the pre-0.8
# <output_dir>/<name><suffix> convention.
output_path <- function(args, id, suffix) {
  if (!is.null(args$outputs[[id]])) return(args$outputs[[id]])
  if (is.null(args$output_dir) || is.null(args$name)) {
    stop(sprintf("no --output %s=PATH and no --output_dir/--name to fall back on", id), call. = FALSE)
  }
  file.path(args$output_dir, paste0(args$name, suffix))
}
