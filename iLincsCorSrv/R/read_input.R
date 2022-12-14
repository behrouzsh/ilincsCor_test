#' @export
read_input <- function(lib, filename) {
  libs<-get("libs", envir=.pkgglobalenv)
  x_lib <- libs[[lib]]
  if (!is.null(x_lib)) {
    return(x_lib$read_input(filename))
  }
}
