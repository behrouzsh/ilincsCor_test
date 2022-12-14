#' Create a iLincsCor Object from the iLincsCor C++ Class
#'
#' Allows for the creation of a iLincsCor Object in _C++_ from _R_
#' using the _C++_ iLincsCor class.
#'
#' @param prefix Path to the library
#'
#' @return
#' A `iLincsCor` object from the _C++_ iLincsCor Class.
#'
#' @examples
#' ##################
#' ## Constructor
#'
#' # Construct new iLincsCor object called "ben"
#' lib <- new(iLincsCor, prefix = "/Users/michal/tmp/data/LIB_5/")
#'
#' ##################
#' ## read input vector
#'
#' v <- lib$read_input("tests/LINCSCP_100.tsv")
#'
#' ## run correlation on 4 CPUs
#'
#' result <- lib$cor(v,4)
#'
#' @name iLincsCor
#' @export iLincsCor
# Load the Rcpp module exposed with RCPP_MODULE( ... ) macro.
loadModule(module = "RcppIlincsCorEx", TRUE)
