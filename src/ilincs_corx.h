#ifndef ILINCS_CORX_H
#define ILINCS_CORX_H

#ifdef RcppExport
// #define DATA_MATRIX_DATATYPE Rcpp::NumericVector
#endif

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>

#include <list>
#include <map>
#include <set>
#include <deque>
#include <stack>
#include <queue>
#include <bitset>
#include <algorithm>
#include <numeric>
#include <utility>
#include <sstream>
#include <span>

#define MEMORY_MAPPED_FILE

#ifndef PRECISION
#define PRECISION double
#endif

// #define LLINT long long int
#define LLINT unsigned long
//===================================================================
//===================================================================
//===================================================================

#define FILENAME "matrix_cp.dat"
#define INFONAME "info.dat"
#define WEIGHTNAME "weight.dat"
#define WEIGHT8NAME "weight_8.dat"
#define DATANAME "data.dat"
#define DATA8NAME "data_8.dat"
#define SIGNAMESFILENAME "signature_names.dat"
#define GENEIDSFILENAME "gene_ids.dat"

#define INPUTEXAMPLEFILENAME "tests/sigFile.tsv"

#define DEBUG(x) x

using namespace std;

#define ABS(a) ((a) < 0 ? - (a) : (a))
#define FOR(i,a,b) for(int i=(a);i<(b);++i)
#define REP(i,n)  FOR(i,0,n)
#define MAX(a,b) (a < b ? b : a)

std::vector<std::vector<std::string> > read_tsv(std::string fname);

#ifndef DATA_MATRIX_DATATYPE

#ifndef MEMORY_MAPPED_FILE
#define DATA_MATRIX_DATATYPE std::vector<PRECISION>
#endif

#ifdef MEMORY_MAPPED_FILE
#define DATA_MATRIX_DATATYPE std::span<PRECISION>
#endif

#endif


typedef struct {
   DATA_MATRIX_DATATYPE weight, data, combined;
   std::vector<PRECISION> weight_sum;
   std::vector<string> signature_names;
   std::vector<int> gene_ids;
   long n_cols, n_rows, precision;
} t_data_matrix;

typedef std::vector<PRECISION> t_input;
typedef std::vector<int> t_input_included;
typedef std::vector<int> t_input_map;
typedef std::vector<PRECISION> t_output;

int read_matrix(std::string prefix, t_data_matrix *data_matrix);

t_output cor(t_input input, t_input input_weights, t_input_included input_included, t_data_matrix *data_matrix,
         long each_n_cols=0, long each_start_col=0);
t_output cor_map(t_input input, t_input input_weights, t_input_map input_map, t_data_matrix *data_matrix,
         long each_n_cols=0, long each_start_col=0);

int read_input(std::vector<int> gene_ids, std::string filename, t_input *input, t_input *input_weights, t_input_included *input_included, int *found_geneids,
               t_input *input_src, t_input *input_weights_src, t_input_map *input_map);

// pthreads
//===================================================================
typedef struct
{
   int n_cols;
   int n_rows;
   int thread_i;
   int number_of_threads;
   t_output *global_results;
   t_input *input;
   t_input *input_weights;
   t_input_map *input_map;
   t_data_matrix *data_matrix;
} t_cor_thread_data;

t_output pcor(t_input input, t_input input_weights, t_input_map input_included, t_data_matrix *data_matrix, int number_of_threads);
t_output pcor_map(t_input input, t_input input_weights, t_input_map input_map, t_data_matrix *data_matrix, int number_of_threads);

std::vector<std::vector<std::string> > read_tsv(std::string fname);
#endif
