/*!
 * \file
 * \brief Program to convert .ens files to .deg files
 * \author Michael Meidlinger
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2016 Michael Meidlinger - All Rights Reserved
 *
 */

#include <itpp/itbase.h>
#include "LDPC_DE.hpp"

using namespace itpp;
using namespace std;

/*!
 *
 * Returns an integer vector containing the indices of linearly independent columns 
 * (relative to the input matrix \c M).
 * \param[in] M           Sparse GF2  matrix with \c ncols columns
 * \param[in] pref_cols   Permutation vector of length \c ncols
 * \return                A vector of length equal to the column rank of \c M 
 *
 * \warning 
 *  This function  does not check wether \c pref_cols is a valid permutation vector (i.e. 
 *  has unique entries within \c 0 and <tt>ncols-1<\tt> )
 * 
 * 
 *
 * */
ivec get_lin_indep_cols(GF2mat_sparse M,  const ivec& pref_cols);


int main(int argc, char **argv){
    // Parse file name
    if(argc!=4){
        cout << "Usage: alist_fullrank infile.alist outfile.alist integer_rand_seed" << endl;
        return EXIT_FAILURE;
    }
    
    string infilename(argv[1]);
    string outfilename(argv[2]);
    string rand_seed_string(argv[3]);
    int rand_seed = stoi(rand_seed_string);
    
    LDPC_Parity H_in(infilename);
        
    GF2mat_sparse H = H_in.get_H();
    GF2mat_sparse Ht = transpose(H);
    
    // randomly select prefered rows to keep
    RNG_reset(rand_seed);
    ivec row_pref = sort_index(randu(H.rows()));
    
    ivec lin_indep_rows = get_lin_indep_cols(Ht, row_pref);
    
    int nvar = H.cols();
    int ncheck = H.rows();
    int ncheck_lin_indep = lin_indep_rows.length();

    cout << "Found " << ncheck_lin_indep << "/" << ncheck 
         <<  " linearly independent rows at original indices " 
         << endl << lin_indep_rows << endl;
   // create new sparse matrix
   GF2mat_sparse Ht_new = GF2mat_sparse(nvar, ncheck_lin_indep);
   
   for(int mm=0; mm < ncheck_lin_indep; mm++){
     Ht_new.set_col(mm, Ht.get_col(lin_indep_rows(mm)));
   }
   
   GF2mat_sparse_alist H_alist;
   H_alist.from_sparse(Ht_new, true); 
   LDPC_Parity H_new(H_alist);
   
   // Get corresponding degree distribution and print it out
   cout << "Calculating LDPC node distributions ..." << endl;
   LDPC_Ensemble ens = get_empirical_ensemble(H_new);
   cout << "Done. Rate of Ensemble = " << ens.get_rate() << endl;
   cout << ens;

   // Save alist
   H_alist.write(outfilename);
   cout << "Successfully saved output file " << outfilename ;

    
   return EXIT_SUCCESS;
}


ivec get_lin_indep_cols(GF2mat_sparse M, const ivec& pref_cols)
{
  int rows = M.rows();
  int cols = M.cols();
  bvec lin_indep_cols = zeros_b(cols); 
  for(int rr=0; rr<rows; rr++){
    // Find cols that contain a '1' in row rr,
    // ordered by preferance
    ivec idx;
    for(int cc=0; cc<cols; cc++){
      if(lin_indep_cols(pref_cols(cc))){ //column already processed and marked as part of set of LI columns
        continue;
      }
      else{
        if( M(rr,pref_cols(cc)) ){
          idx = concat(idx, pref_cols(cc));
        }
      }
     
    }   
    if( length(idx)>0 ){
      //Choose the row according to preferance
      int kk = idx(0);
      lin_indep_cols(kk) = true;
      idx.del(0);
      // Subtract col from other cols
      for(int ii=0; ii<length(idx); ii++){
            M.set_col(idx(ii),  M.get_col(idx(ii)) + M.get_col(kk) );      
      }
    }
  }
  return find(lin_indep_cols);

}

