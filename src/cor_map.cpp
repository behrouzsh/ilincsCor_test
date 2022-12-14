
#ifdef ILINCS_RCPP
#include <Rcpp.h>
#define DATA_MATRIX_DATATYPE Rcpp::NumericVector
#endif

#include "ilincs_corx.h"

t_output cor_map(t_input input, t_input input_weights, t_input_map input_map, t_data_matrix *data_matrix,
         long each_n_cols, long each_start_col)
{
   // cout << "cor: matrix: weights " << data_matrix->weight[0] << ", " << data_matrix->weight[1]  << endl;

   struct timeval start_time;
   gettimeofday(&start_time, NULL);

   t_output each_result;
   long n_rows=input_map.size();
   long n_cols=data_matrix->n_cols;
   if (n_cols==0 || n_rows==0) 
   {
      // cerr << "Matrix not loaded" << endl;
      return each_result;
   }

   long each_start=each_start_col*n_rows;

   if (each_n_cols == 0) 
   {
      each_n_cols=n_cols;
   }

   each_result.resize(each_n_cols);

   std::vector<PRECISION> weight_cache(n_rows);

   LLINT block = each_start;

   for (long i = 0; i < each_n_cols;  i++)
   {
      // compute coefs
      PRECISION sum_wixi = 0.0;
      PRECISION sum_wiyi = 0.0;
      PRECISION sum_wi   = 0.0;

      // TODO: if weights are not provided the extra addition is eliminted and speed up happens
      // no user weights:: CAN'T HAVE precomputed variables cov_xxw_numerator, sum_wi, sum_wixi
      // no user weights:: since we don't know which elements are included

      // omp on this level makes things lot worse
      for (long j = 0; j < n_rows; j++) 
      {
	 PRECISION _libvec = data_matrix->data[block + input_map[j]]; // xi
         PRECISION _input = input[j];                                 // yi
         PRECISION _weight = (data_matrix->weight[block + input_map[j]] + input_weights[j]); // wi
         // PRECISION _weight = input_included[j] * (data_matrix->weight[block + j]); // no weights

         sum_wixi += _weight * _libvec;
         sum_wiyi += _weight * _input; // if not included _input = 0
         sum_wi   += _weight;

	 // record weight so we don't have to check input_included or ->weight in the next for-loop
	 weight_cache[j] = _weight;
      }

      PRECISION m_xw = sum_wixi / sum_wi;
      PRECISION m_yw = sum_wiyi / sum_wi;

      PRECISION cov_xxw_numerator = 0.0;
      PRECISION cov_xyw_numerator = 0.0;
      PRECISION cov_yyw_numerator = 0.0;

      // omp on this level makes things lot worse
      for (long j = 0;  j < n_rows;  j++)
      {
         PRECISION _libvec  = ( data_matrix->data[block + input_map[j]] - m_xw );
         PRECISION _input   = ( input[j]                                - m_yw );
         PRECISION _weight  = weight_cache[j];

         cov_xxw_numerator += (_libvec ) * (_libvec ) * (_weight );
         cov_xyw_numerator += (_libvec ) * (_input  ) * (_weight );
         cov_yyw_numerator += (_input  ) * (_input  ) * (_weight );
      }

      each_result[i] = cov_xyw_numerator / sqrt(cov_yyw_numerator * cov_xxw_numerator);

      block += n_rows;
   }

   struct timeval end_time;
   double delta_time;
   gettimeofday(&end_time, NULL);
   delta_time = (double) (end_time.tv_sec - start_time.tv_sec);
   delta_time += (double) (end_time.tv_usec - start_time.tv_usec) / 1.0e+6;

   // cerr << "cor finished (" << delta_time << "s)" << endl;
   // cerr << "output: " << each_result[0] << ", " << each_result[1] << endl;
   return(each_result);
}

