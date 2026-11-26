# Phase timings as obkit structured events (<dir>/obkit-events.jsonl).
#
# timing_init(events_path) points obkit at the directory of the declared
# events output; obkit always names the file obkit-events.jsonl.
# with_phase("pca", expr) emits a start event, evaluates expr, emits the end
# event and returns expr's value. `attrs` (named list, or a function of the
# value returning one) goes on the end event, e.g. the number of clustering runs.

timing_init <- function(events_path) {
  if (basename(events_path) != "obkit-events.jsonl") {
    stop("the events output must be named obkit-events.jsonl, got ", events_path, call. = FALSE)
  }
  obkit::logger_init(dirname(events_path))
}

with_phase <- function(name, expr, attrs = NULL) {
  obkit::logger_emit(name, "start")
  value <- expr
  if (is.function(attrs)) attrs <- attrs(value)
  obkit::logger_emit(name, "end", attrs)
  invisible(value)
}
