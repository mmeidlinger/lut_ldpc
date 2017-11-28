/*!
 * \file
 * \brief Program to generate LDPC Codes
 * \author Michael Meidlinger
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2016 Michael Meidlinger - All Rights Reserved
 *
 */

#include "LDPC_DE.hpp"
#include "LDPC_Code_LUT.hpp"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <boost/format.hpp>

using namespace itpp;
using namespace std;
int main(int argc, char **argv){

    // Generate Parity and optionally, a corresponding Generator
    
    bool make_generator = false;
    bool make_bp_codec = false; //specify parameters below
    bool make_lut_codec = true; //specify parameters below
    bool save_parity = true;
    
    
    

    // Set Codeword size and edge distributions
    int N = 1e5;
    
    // Girth target.
    int girth_target = 0;
    
    // Check for large code sizes when using a generator
    if(make_generator && N>=1e5){
        it_error("For Codelengths larger than 1e5, a generator matrix uses too much memory!");
    }
    
    // Set optimization based on codeword length
    std::string opt_goals;
    if(N<=0)
        it_error("Codeword length must be positive!");
    else if(N<=5000)
        opt_goals = "500 16";
    else if(N<=50000)
        opt_goals = "150 8";
    else if(N<=500000)
        opt_goals = "100 16";
    else
        opt_goals = "0 0";
        

            
    
    // ensemble filename
   // std::string ensemble_filename = "ensembles/rate0.50_dv02-09_dc06-08.ens";
    std::string ensemble_filename = "ensembles/rate0.50_dv02-17_dc08-09_lut_q4.ens";
  //  std::string ensemble_filename = "ensembles/rate0.50_dv03_dc06.ens";
    
    // Code filename without suffix
   // std::string code_filename = "codes/rate0.50_irreg_dv02-09_dc06-08";
     std::string code_filename = "codes/rate0.50_dv02-17_dc08-09_lut_q4";
   // std::string code_filename = "codes/rate0.50_dv03_dc06";
    code_filename = code_filename + "_N" + to_str(N);
    
    
    LDPC_Ensemble ens(ensemble_filename);
    
    cout << "Designing code of length " << N << endl;
    cout << "Design Rate: " << ens.get_rate() << endl;
    cout << ens.get_var_degree_dist() << endl << ens.get_chk_degree_dist() << endl;
    cout << "Optimization Goals: " << opt_goals << endl
         << "Girth target: " << girth_target << endl;
    

    it_info("Creating Parity...");
    LDPC_Parity_Irregular H(N, ens.get_var_degree_dist(),
                            ens.get_chk_degree_dist(),
                            "rand",
                            opt_goals);
    
    it_info("Girth optimization with target " << girth_target);
    H.cycle_removal_MGW(girth_target);
    it_info("done.");
    
    if(save_parity)
        H.save_alist( code_filename + ".alist" );
    
    LDPC_Generator_Systematic G;
    if(make_generator){
        it_info("Making a Generator Matrix ...");
        G.construct(&H);
        G.save( code_filename + ".gen.it");
    }
    
    if(make_bp_codec){
        it_info("Building BP Codec ...");
        LDPC_Code C;
        if(make_generator){
            C = LDPC_Code(&H, &G);
        }
        else{
            C = LDPC_Code(&H);
        }
        int max_iters = 200;
        C.set_exit_conditions(max_iters);
        C.save_code(code_filename + ".bp_codec.it");
    }
    
    if(make_lut_codec){
        it_info("Building LUT Codec ...");
        // Set parameters
        std::string tree_method = "auto_bin_balanced"; // auto_bin_balanced, auto_bin_high, filename=<filename>
        LDPC_Ensemble ens = get_empirical_ensemble(H);
        cout << "Empirical Rate and ensemble: " << ens.get_rate() << endl;
        cout << ens.get_var_degree_dist() << endl << ens.get_chk_degree_dist() << endl;
        
        int max_iter = 100;
        bool min_lut = true;
        ivec Nq_Msg = ones_i(max_iter)* pow2i(4);
        int Nq_Cha = pow2i(4);
        bvec reuse_vec = zeros_b(max_iter);
        double design_thr = 0.864395;
        LDPC_Code_LUT C;
        if(make_generator){
            C = LDPC_Code_LUT(&H, &G);
        }
        else{
            C = LDPC_Code_LUT(&H);
        }
        
        design_thr = C.design_luts(tree_method, ens, min_lut, sqr(design_thr), max_iter, reuse_vec, Nq_Cha, Nq_Msg);
        C.save_code(code_filename + "_" + tree_method +
                    "_MinLUT" + to_str(min_lut) +
                    "_maxIter" + to_str(max_iter) +
                    "_designThr" + boost::str(boost::format("%3g") % design_thr) +
                    "_NqCha" + to_str(Nq_Cha) +
                    "_NqMsg" + to_str(itpp::min(Nq_Msg)) + "-" + to_str(itpp::max(Nq_Msg)) +
                    ".lut_codec.it");
    }
    
    

    it_info("Code Generation completed successfully.");
    



    return EXIT_SUCCESS;
}
