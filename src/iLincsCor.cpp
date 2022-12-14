#include <Rcpp.h>
#include <iostream>
using namespace std;

#include "iLincsCor.h"

using namespace Rcpp;

// $read_matrix(prefix)
// $read_input(filename)
// $cor(input,cpus)

// Constructor
iLincsCor::iLincsCor(std::string prefix) {
   cout << "constructor iLincsCor" << endl;
   int res = read_matrix(prefix, &(this->data_matrix));
   this->genes = this->data_matrix.n_rows;
   this->signatures = this->data_matrix.n_cols;
   this->signature_names = this->data_matrix.signature_names;
   this->gene_ids = this->data_matrix.gene_ids;
}

List iLincsCor::R_read_input(std::string input_filename){

  cout << "R_read_input(" << input_filename << ")" << endl;

  int found_geneids=0;

  // preallocate arrays to match the library
  cout << "Preallocate:" << data_matrix.n_rows << endl;
  t_input input(data_matrix.n_rows);
  t_input input_weights(data_matrix.n_rows);
  t_input_included input_included(data_matrix.n_rows);

  t_input input_src(0);
  t_input input_weights_src(0);
  t_input_included input_map(0);

  int res = read_input(data_matrix.gene_ids, input_filename, &input, &input_weights, &input_included, &found_geneids, &input_src, &input_weights_src, &input_map);
  if (res == 1 ) {
     return 1;
  }

  // cout << "size=" << vec_size << endl;
  // cout << "path=" << path << endl;
  // List df = List::create ();

  List df = List::create(
        Named("found_geneids") = found_geneids,
        _("input") = input,
        _("input_weights") = input_weights,
        _("input_included") = input_included,
        _("input_src") = input_src,
        _("input_weights_src") = input_weights_src,
        _("input_map") = input_map
        );

  return df;
}

Rcpp::NumericVector iLincsCor::R_cor(List input, int workers){
  // cout << "R_cor(" << workers << ")" << endl;

   t_input input_data = as<t_input >(input["input"]);
   t_input input_weights = as<t_input >(input["input_weights"]);
   t_input_included input_included = as<t_input_included >(input["input_included"]);

   // cout << "matrix: " << data_matrix.n_cols << " x " << data_matrix.n_rows << endl;
   // cout << "weights front: " << data_matrix.weight[0] << endl;

   t_output output;
   // output = cor(input, input_weights, input_included, &data_matrix, 0, 0);
   output = pcor(input_data, input_weights, input_included, &data_matrix, workers);

   return Rcpp::wrap(output);
}

Rcpp::NumericVector iLincsCor::R_cor_map(List input, int workers){
  // cout << "R_cor(" << workers << ")" << endl;

   t_input input_src = as<t_input >(input["input_src"]);
   t_input input_weights_src = as<t_input >(input["input_weights_src"]);
   t_input_map input_map = as<t_input_map >(input["input_map"]);

   // cout << "matrix: " << data_matrix.n_cols << " x " << data_matrix.n_rows << endl;
   // cout << "weights front: " << data_matrix.weight[0] << endl;

   t_output output;

   output = pcor_map(input_src, input_weights_src, input_map, &data_matrix, workers);

   return Rcpp::wrap(output);
}

