# Load only the function definitions from analysis.R (no data loading, no plots),
# so the tests exercise the real code without running the whole script.
suppressPackageStartupMessages(library(ggplot2))

load_analysis_functions <- function(env = parent.frame()) {
  path <- file.path(project_root(), "analysis.R")
  exprs <- parse(path, keep.source = FALSE)
  for (e in exprs) {
    is_fn_def <- is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))
    if (is_fn_def) eval(e, envir = env)
  }
  invisible(env)
}

project_root <- function() {
  p <- normalizePath(getwd())
  while (!file.exists(file.path(p, "analysis.R")) && dirname(p) != p) p <- dirname(p)
  p
}

load_analysis_functions(globalenv())
