#' @export
signature_names <- function(lib) {
  libs<-get("libs", envir=.pkgglobalenv)
  x_lib <- libs[[lib]]
  if (!is.null(x_lib)) {
    return(x_lib$signature_names)
  }
}
