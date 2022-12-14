#' @import iLincsCor
.onLoad <- function(libname, pkgname){
  library(iLincsCor)
  x_datadir <- Sys.getenv("ILINCSCOR_DATADIR")
  # if (x_datadir == "") x_datadir <- "/opt/raid10/docker/ilincs_cor/data"
  if (x_datadir == "") x_datadir <- "/data"

  x_names <- list.dirs(path = x_datadir, full.names = FALSE, recursive = FALSE)
  x_dirs <- paste0(list.dirs(path = x_datadir, full.names = TRUE, recursive = FALSE), "/")
  x_libs <- sapply(x_dirs, function(x_dir) { new(iLincsCor,x_dir) })
  names(x_libs) <- x_names
  assign(x="libs", value=x_libs, envir=.pkgglobalenv)

  # assign(url, file, envir=cacheEnv)
  # get(url, envir=cacheEnv)
}

.pkgglobalenv <- new.env(parent=emptyenv())
