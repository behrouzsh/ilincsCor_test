#include "ilincs_corx.h"
#include "util.h"

inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

/*
 * Read tsv input file and for a given library format the vector of value, weight and whether the gene is included
 * also return the number of gene ids found
 */
int read_input(std::vector<int> gene_ids, std::string filename, t_input *input, t_input *input_weights, t_input_included *input_included, int *found_geneids,
               t_input *input_src, t_input *input_weights_src, t_input_map *input_map)
{
   *found_geneids=0;

   // Read input
   std::vector<std::vector<std::string> > input_vec_vec;
   input_vec_vec = read_tsv(filename.c_str());
   int n = input_vec_vec.size();

   if (n == 0) return -1;

   int geneid_column=-1; // default?
   int LogDiffExp_column=-1; // default?
   int pvalue_column=-1; // default?
   int header_size = input_vec_vec[0].size();

   // find which column is geneid, logdiffexp and pvalue (if exists)
   for(int i = 0; i < header_size; i++) {
         DEBUG(cerr << "Checking column " << input_vec_vec[0][i] << endl;)
         if (ends_with(input_vec_vec[0][i],"_geneid")) {
            geneid_column=i;
            DEBUG(cerr << "Found geneid at " << i << endl;)
         }
         if (ends_with(input_vec_vec[0][i],"_LogDiffExp")) {
            LogDiffExp_column=i;
            DEBUG(cerr << "Found LogDiffExp at " << i << endl;)
         }
         if (ends_with(input_vec_vec[0][i],"_pvalue")) {
            pvalue_column=i;
            DEBUG(cerr << "Found pvalue at " << i << endl;)
         }
   }

   // we must have logdiffexp
   if (LogDiffExp_column==-1) {
      DEBUG(cerr << "logdiffexp column not found" << endl;)
      return -1; 
   }

   // read input vector and assign logdiffexp to appropriate geneids positions
   std::fill((*input).begin(), (*input).end(), 0);
   std::fill((*input_weights).begin(), (*input_weights).end(), 0);
   std::fill((*input_included).begin(), (*input_included).end(), 0);
   int n_rows = input->size();
   for(int i = 1; i < n && i <= n_rows; i++){
      // DEBUG(cerr << "Converting " << input_vec_vec[i][LogDiffExp_column] << " at line " << i << " column " << LogDiffExp_column << endl;)

      int gene_position=i-1;

      PRECISION input_exp;
      try {
          input_exp = std::stod(input_vec_vec[i][LogDiffExp_column]);
      }
      catch (int e) {
         DEBUG(cerr << "Can't convert diffexp " << input_vec_vec[i][LogDiffExp_column] << " at line " << i << " column " << LogDiffExp_column << endl;)
         return -1;
      }

      PRECISION input_pvalue=0;
      if (pvalue_column!=-1) {
         try {
            input_pvalue = std::stod(input_vec_vec[i][pvalue_column]);
            input_pvalue = -1*log10(input_pvalue);
         }
         catch (int e) {
            DEBUG(cerr << "Can't convert pvalue " << input_vec_vec[i][pvalue_column] << " at line " << i << " column " << pvalue_column << endl;)
            return -1;
         }
      }
     
      if (geneid_column!=-1) {
         int gene_id;
         try {
            gene_id = std::stod(input_vec_vec[i][geneid_column]);
         } 
         catch (int e) {
            DEBUG(cerr << "Can't convert gene_id " << input_vec_vec[i][geneid_column] << " at line " << i << " column " << geneid_column << endl;)
            return -1;
         }

         // shortcut
         if (gene_id != gene_ids[gene_position]) {

            // FIXME: make this binary search
            std::vector<int>::iterator it = std::find(gene_ids.begin(), gene_ids.end(), gene_id);
            if (it == gene_ids.end()) {
               gene_position = -1; // can't find gene_id
            } else {
               gene_position = std::distance(gene_ids.begin(), it);
            }
         }
      }
      if (gene_position != -1) {
         if ((*input_included)[gene_position]==0) {
            (*found_geneids)++;
         }
         (*input)[gene_position]=input_exp;
         (*input_weights)[gene_position]=input_pvalue;
         (*input_included)[gene_position]=1;

	 // map
	 if (input_src) input_src->push_back(input_exp);
	 if (input_weights_src) input_weights_src->push_back(input_pvalue);
	 if (input_map) input_map->push_back(gene_position);
      }
   }
   return 0;
}
