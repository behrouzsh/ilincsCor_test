#ifndef ILINCS_COR_H
#define ILINCS_COR_H

#include <string>
#include <vector>

#include <Rcpp.h>
using namespace Rcpp;

#include "ilincs_corx.h"

class iLincsCor
{
public:
    t_data_matrix data_matrix; // Must be hidden from R
    int genes;
    int signatures;
    std::vector<string> signature_names;
    std::vector<int> gene_ids;

    // Constructor
    iLincsCor(std::string prefix);

    // Exposed functions
    // int R_read_matrix(std::string prefix);
    List R_read_input(std::string input_filename);
    Rcpp::NumericVector R_cor(List input, int workers);
    Rcpp::NumericVector R_cor_map(List input, int workers);
};

#endif /* ILINCS_COR_H */
