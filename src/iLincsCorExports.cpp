#include <Rcpp.h>

#include "iLincsCor.h"

RCPP_MODULE(RcppIlincsCorEx) {
    using namespace Rcpp;
    class_<iLincsCor>("iLincsCor")
        .constructor<std::string>()
      .method("read_input", &iLincsCor::R_read_input)
      .method("cor", &iLincsCor::R_cor)
      .method("cor_map", &iLincsCor::R_cor_map)

      .field("genes", &iLincsCor::genes)
      .field("signatures", &iLincsCor::signatures)
      .field("gene_ids", &iLincsCor::gene_ids)
      .field("signature_names", &iLincsCor::signature_names)
    ;
}

