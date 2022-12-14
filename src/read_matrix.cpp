#ifdef ILINCS_RCPP
#include <Rcpp.h>
#define PRECISION double
// we use std::vector<double> and assume implicit converstion with Rcpp
// for MEMORY_MAPPED_FILE we use std::span
// #define DATA_MATRIX_DATATYPE Rcpp::NumericVector
#endif

#include "ilincs_corx.h"

int read_matrix(std::string prefix, t_data_matrix *data_matrix)
{
   struct timeval start_time;
   gettimeofday(&start_time, NULL);

   std::string info_filename = prefix + INFONAME;
   std::string weight_filename = prefix + WEIGHTNAME;
   std::string data_filename = prefix + DATANAME;
   std::string signature_names_filename = prefix + SIGNAMESFILENAME;
   std::string gene_ids_filename = prefix + GENEIDSFILENAME;

   if (sizeof(PRECISION)==sizeof(double)) {
      weight_filename = prefix + WEIGHT8NAME;
      data_filename = prefix + DATA8NAME;
   }

   // Reading matrix ======================================================================
   int fd_weight = open(weight_filename.c_str(), O_RDONLY);
   int fd_data = open(data_filename.c_str(), O_RDONLY);
   ifstream ifs_info(info_filename.c_str());
   ifstream ifs_signature_names(signature_names_filename.c_str());
   ifstream ifs_gene_ids(gene_ids_filename.c_str());

   if (fd_weight==-1 || fd_data==-1 || ifs_info.fail()) {
      return 1;
   }

   // read cols/rows from info
   std::string line;
   while(getline(ifs_info, line)){
      line.erase(std::remove_if(line.begin(), line.end(), ::isspace),
            line.end());
      if(line[0] == '#' || line.empty())
         continue;
      int delimiterPos = line.find("=");
      std::string name = line.substr(0, delimiterPos);
      std::string value = line.substr(delimiterPos + 1);
      DEBUG(std::cerr << name << " ... " << value << '\n';)
      if (name.compare("signatures")==0) data_matrix->n_cols=stoi(value);
      if (name.compare("genes")==0) data_matrix->n_rows=stoi(value);
      // if (name.compare("precision")==0) data_matrix->precision=stoi(value);
   }

   if (data_matrix->n_cols==0 || data_matrix->n_rows==0) {
      DEBUG(cerr << "Error: info file missing data (" << data_matrix->n_cols << "," << data_matrix->n_rows << ")" << endl;)
      return 1;
   }

   LLINT NUM_ELEMENTS = data_matrix->n_rows*data_matrix->n_cols;

   DEBUG(cerr << "Matrix: " << data_matrix->n_rows << " x " << data_matrix->n_cols << endl;)
   DEBUG(cerr << "Precision: " << sizeof(PRECISION) << endl;)
   DEBUG(cerr << "Memory (2*" << NUM_ELEMENTS*sizeof(PRECISION)/1024/1024 << "MB)" << endl;)
   // data_matrix->weight.resize(NUM_ELEMENTS);
   // data_matrix->data.resize(NUM_ELEMENTS);


#pragma omp parallel sections
   {
#pragma omp section
      {
         if (!ifs_signature_names.fail()) {
            string line;
            while(getline(ifs_signature_names, line)) {
               // remove quotes
               if ( line.front() == '"' ) {
                  line.erase( 0, 1 ); // erase the first character
                  line.erase( line.size() - 1 ); // erase the last character
               }
               data_matrix->signature_names.push_back(line);
            }
         }
         // cout << "read_matrix: signature_names loaded... " << endl;
      }

#pragma omp section
      {
         if (!ifs_gene_ids.fail()) {
            string line;
            while(getline(ifs_gene_ids, line)) {
               // remove quotes
               if ( line.front() == '"' ) {
                  line.erase( 0, 1 ); // erase the first character
                  line.erase( line.size() - 1 ); // erase the last character
               }
               data_matrix->gene_ids.push_back(std::stoi(line));
            }
            // create invert of gene_ids for fast search
         }
         // cout << "read_matrix: gene_ids loaded... " << endl;
      }

#pragma omp section
      {
#ifdef MEMORY_MAPPED_FILE
// #warning "Memory Mapped File"
         // cout << "read_matrix: data setting up memory mapped file... " << endl;
	 PRECISION* b = (PRECISION*)(mmap(NULL, sizeof(PRECISION)*NUM_ELEMENTS, PROT_READ, MAP_PRIVATE, fd_weight, 0u));
         data_matrix->weight={b,NUM_ELEMENTS};
         DEBUG(cerr << "-- weight " << data_matrix->weight.size() << endl;)
#endif

#ifndef MEMORY_MAPPED_FILE
   	 data_matrix->weight.resize(NUM_ELEMENTS);
         LLINT mem_counter=0;
         LLINT mem_size=NUM_ELEMENTS;
         LLINT max_chunk=2147483647/sizeof(PRECISION);
         while(mem_counter<mem_size) {
            LLINT chunk=(mem_size-mem_counter);
            if (chunk>max_chunk) chunk=max_chunk;
            // cerr << "Reading weight ... " << mem_counter << "/" << mem_size << " (" << chunk << ")" << endl;
            read(fd_weight, &(data_matrix->weight[mem_counter]), chunk*sizeof(PRECISION));
            mem_counter+=chunk;
         }
#endif

         // cerr << "weights: " << data_matrix->weight[0] << ", " << data_matrix->weight[1] << "... " << endl;
      }
#pragma omp section
      {
#ifdef MEMORY_MAPPED_FILE
	 PRECISION* b = (PRECISION*)(mmap(NULL, sizeof(PRECISION)*NUM_ELEMENTS, PROT_READ, MAP_PRIVATE, fd_data, 0u));
         data_matrix->data={b,NUM_ELEMENTS};
#endif

#ifndef MEMORY_MAPPED_FILE
         data_matrix->data.resize(NUM_ELEMENTS);
         LLINT mem_counter=0;
         LLINT mem_size=NUM_ELEMENTS;
         LLINT max_chunk=2147483647/sizeof(PRECISION);
         while(mem_counter<mem_size) {
            LLINT chunk=(mem_size-mem_counter);
            if (chunk>max_chunk) chunk=max_chunk;
            // cerr << "Reading data ... " << mem_counter << "/" << mem_size << " (" << chunk << ")" << endl;
            read(fd_data, &(data_matrix->data[mem_counter]), chunk*sizeof(PRECISION));
            mem_counter+=chunk;
         }
#endif

         // cerr << "data: " << data_matrix->data[0] << ", " << data_matrix->data[1] << "... " << endl;
      }
   }
/*
   if (0)
   {
     data_matrix->weight_sum.resize(data_matrix->n_cols);
     for(LLINT i=0;i<data_matrix->n_cols;i++) {
       data_matrix->weight_sum[i] = std::accumulate(data_matrix->weight.begin()+i*data_matrix->n_rows, data_matrix->weight.begin()+(i+1)*data_matrix->n_rows-1, 0);
     }
   }

   if (0)
   {
         data_matrix->combined.resize(2*NUM_ELEMENTS);
	 PRECISION *combined = &(data_matrix->combined[0]);
	 PRECISION *data = &(data_matrix->data[0]);
	 PRECISION *weight = &(data_matrix->weight[0]);

         for(LLINT i=0;i<NUM_ELEMENTS;i++) {
	    combined[i*2]=data[i];
	    combined[i*2+1]=weight[i];
	 }
         cerr << "combined: " << data_matrix->combined[0] << ", " << data_matrix->combined[1] << "... " << endl;
   }
*/
   struct timeval end_time;
   double delta_time;
   gettimeofday(&end_time, NULL);
   delta_time = (double) (end_time.tv_sec - start_time.tv_sec);
   delta_time += (double) (end_time.tv_usec - start_time.tv_usec) / 1.0e+6;

   // cerr << "Return read_matrix (" << delta_time << "s)" << endl;
   return 0;
}
