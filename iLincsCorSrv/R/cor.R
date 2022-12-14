#' @export
cor <- function(lib, sig, n=4) {
  libs<-get("libs", envir=.pkgglobalenv)
  x_lib <- libs[[lib]]
  if (!is.null(x_lib)) {
    return(x_lib$cor(sig,n))
  }
}
