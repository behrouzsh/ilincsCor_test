#include "ilincs_corx.h"

void *cor_thread_fn(void *cor_thread_data);

t_output pcor(t_input input, t_input input_weights, t_input_included input_included, t_data_matrix *data_matrix, int number_of_threads) 
{
   struct timeval start_time;
   gettimeofday(&start_time, NULL);

   pthread_t *running_thread = (pthread_t*)malloc(number_of_threads*sizeof(pthread_t));
   t_cor_thread_data *cor_thread_params = (t_cor_thread_data*)malloc(number_of_threads*sizeof(t_cor_thread_data));

   // Create threads ================================================================
   t_output results(data_matrix->n_cols);

   for (int thread_i = 0; thread_i < number_of_threads; thread_i++)
   {
      cor_thread_params[thread_i].n_cols = data_matrix->n_cols;
      cor_thread_params[thread_i].n_rows = data_matrix->n_rows;
      cor_thread_params[thread_i].thread_i = thread_i;
      cor_thread_params[thread_i].number_of_threads = number_of_threads;
      cor_thread_params[thread_i].global_results = &results;
      cor_thread_params[thread_i].input = &input;
      cor_thread_params[thread_i].input_weights = &input_weights;
      cor_thread_params[thread_i].input_map = &input_included;
      cor_thread_params[thread_i].data_matrix = data_matrix;

      pthread_create(&running_thread[thread_i], NULL, cor_thread_fn, &cor_thread_params[thread_i]);
   }

   for (int thread_i = 0; thread_i < number_of_threads; thread_i++)
   {
        pthread_join(running_thread[thread_i], NULL);
   }

   free(cor_thread_params);
   free(running_thread);

/*
   struct timeval end_time;
   double delta_time;
   gettimeofday(&end_time, NULL);
   delta_time = (double) (end_time.tv_sec - start_time.tv_sec);
   delta_time += (double) (end_time.tv_usec - start_time.tv_usec) / 1.0e+6;

   cerr << "pcor finished (" << delta_time << "s)" << endl;
*/
   return results;
}

void *cor_thread_fn(void *cor_thread_data) {
   int n_cols=((t_cor_thread_data*)cor_thread_data)->n_cols;
   // int n_rows=((t_cor_thread_data*)cor_thread_data)->n_rows;
   int thread_i=((t_cor_thread_data*)cor_thread_data)->thread_i;
   int num_threads=((t_cor_thread_data*)cor_thread_data)->number_of_threads;
   t_output *global_results=((t_cor_thread_data*)cor_thread_data)->global_results;
   t_input *input=((t_cor_thread_data*)cor_thread_data)->input;
   t_input *input_weights=((t_cor_thread_data*)cor_thread_data)->input_weights;
   t_input_map *input_map=((t_cor_thread_data*)cor_thread_data)->input_map;
   t_data_matrix *data_matrix = ((t_cor_thread_data*)cor_thread_data)->data_matrix;

   int each_start_col = (thread_i * n_cols) / num_threads;
   int each_start_col_n1 = ((thread_i+1) * n_cols) / num_threads;
   int each_n_cols = each_start_col_n1-each_start_col;

   t_output result = cor(*input, *input_weights, *input_map, data_matrix, each_n_cols, each_start_col);

   for ( int i = 0;  i < each_n_cols;  i++)
   {
      (*global_results)[each_start_col + i] = result[i];
   }
   return NULL;
}

