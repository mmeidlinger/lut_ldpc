/*!
 * \file
 * \brief Program to optimize degree distributions
 * \author Michael Meidlinger
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2016 Michael Meidlinger - All Rights Reserved
 *
 */

#include "LDPC_DE.hpp"
#include <boost/program_options.hpp>
#include <itpp/itbase.h>
#include <thread>


#define MAX_LLR_MAGNITUDE 25
#define MAX_BISEC_ITER 50
#define NQ_FINE 5000
#define PE_MAX 1e-17
#define THR_PREC 1e-7

namespace po = boost::program_options;
using namespace lut_ldpc;
using namespace std;
using namespace itpp;

/*!
 \brief Wrapper to evolve an ensemble with noise level \c thr for multithreaded execution
 */
void evolve_thread(LDPC_DE_LUT de, double thr, double Pe_max, double* Pe, int* iters){
    mat Pmat_dummy;
    vec Pe_trace;
    de.evolve(thr, true, false, Pmat_dummy, Pe_trace);
    ivec pe_idx_vec = find(Pe_trace < Pe_max); 
    *Pe = Pe_trace.right(1).get(0);

    if(pe_idx_vec.length() > 0){
        *iters = pe_idx_vec.get(0);
    }
}



int main(int argc, char **argv){

    cout << "Called program via " << endl << argv[0];
    for(int ii=1; ii<argc; ii++){
        cout << " " << argv[ii];
    } 
    cout << endl;

    bool min_approx = false;   // true if min-LUT algorithm should be used
    int Nq_msg, Nq_cha;
    int maxiter;
    double thr;
    string ensemble_filename;
    double scale_down;
    LDPC_Ensemble ens;
    std::vector<string> degree_dist;
    std::vector<int> reuse_vec_in;
    double Pe_max;
    int max_reuse_stages;
    string lut_tree_design, lut_table_design;

    try {
        po::options_description desc(
                    "OPTIONS:");
        desc.add_options()
        ("min-approx,m",
            po::bool_switch(&min_approx),
	 		"Approximate Check Node Updates" )
        ("quant-bits-msg",
            po::value<int>(&Nq_msg)->default_value(4),
            "Number of quantization bits for messages")
        ("quant-bits-cha",
            po::value<int>(&Nq_cha)->default_value(4),
            "Number of quantization bits for channel outputs")
        ("threshold,t",
            po::value<double>(&thr)->required(),
            "Noise value to run DE. If not provided, found by bisction bisec_search")
        ("ensemble,e",
            po::value<string>(&ensemble_filename),
            "Filename for initial ensemble")
        ("iterations,i",
            po::value<int>(&maxiter)->default_value(100),
            "Number of Message passing iterations")
        ("degree-dist,d",
            po::value<std::vector<string> >(&degree_dist)->multitoken(),
            "Degree ditribution is the form,  \"VN_degrees / VN_probabilities / CN_degrees / CN_probabilities \"")
        ("scale-down,s",
            po::value<double>(&scale_down)->default_value(0.995),
            "Scale down threshold by this value if an updated ensemble fails to converge")
        ("pmax,p",
            po::value<double>(&Pe_max)->default_value(1e-11),
            "Convergence error probability")
        ("reuse-stages,r",
            po::value<int>(&max_reuse_stages)->required(),
            "Number of distinct LUT stages")
        ("reuse-vec,v",
            po::value<std::vector<int>>(&reuse_vec_in)->multitoken(),
            "Provide an initial reuse vector")
        ("lut-table-design",
            po::value<string>(&lut_table_design)->default_value("joint_root"),
            "Strategy for LUT table creation")
        ("lut-tree-design",
            po::value<string>(&lut_tree_design)->default_value("auto_bin_balanced"),
            "Strategy for LUT table creation")
        ("help,h", "produce help message")
        ;
        

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        
        if (vm.count("help")) {
            cout << desc << "\n";
            return EXIT_SUCCESS;
        }
        
        // We put this after processing the help message to avoid required arguments if all we need is "help"
        po::notify(vm);

        //==== Parse degree distribution
        if (vm.count("degree-dist")==0 && vm.count("ensemble")==1) {
            cout << "Reading initial enemble from file" << ensemble_filename << " ... ";
            ens = LDPC_Ensemble(ensemble_filename);
            cout << "Done." << endl;
        }
        else if(vm.count("degree-dist")==1 && vm.count("ensemble")==0){
            cout << "Reading initial enemble from command line ... ";
            
            string delimiter = "/";
            size_t pos = 0;
            int ii = 0;
            string degree_dist_string = "";
            for(unsigned int ii=0; ii<degree_dist.size(); ii++){
                degree_dist_string = degree_dist_string + " " + degree_dist[ii];
            }
            Array<string> degree_dist_array(4);
            while ((pos = degree_dist_string.find(delimiter)) ) {
                degree_dist_array(ii) = degree_dist_string.substr(0, pos);
                degree_dist_string.erase(0, pos + delimiter.length());
                ii++;
                if(pos == std::string::npos) {break;}
                it_assert(ii<=4, "Error reading degree distribution from command line!");
            }
            ivec dv = degree_dist_array(0);
            vec lam = degree_dist_array(1);
            ivec dc = degree_dist_array(2);
            vec rho = degree_dist_array(3);

            ens = LDPC_Ensemble(dv, lam, dc, rho);
            cout << "Done." << endl;

        }
        else{
            it_error("Either --degree-dist or --enseble is required one time");
        }
        // Echo  ensemble
        cout << "Successfully read ensemble with rate " <<  ens.get_rate() << endl << ens << endl;

    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

   


    int Nq_Cha = pow2i(Nq_cha);
    ivec Nq_Msg_vec =   pow2i(Nq_msg)*ones_i(maxiter);
    bvec reuse_vec;
    Array<Array<LUT_Tree>> var_luts;
    Array<Array<LUT_Tree>> chk_luts; 
    get_lut_tree_templates(lut_tree_design, ens, Nq_Msg_vec, Nq_Cha, min_approx, var_luts, chk_luts);
    // Get reuse vector
    if(reuse_vec_in.size() == 0){
        reuse_vec = zeros_b(maxiter);
    }
    else if((int)reuse_vec_in.size() == maxiter){
        reuse_vec.set_size(maxiter);
        for(int ii=0; ii<maxiter; ii++){
            reuse_vec(ii) = reuse_vec_in[ii];
        }
        cout << "Provided initial reuse_vec = " << reuse_vec << endl;
    }
    else{
        it_error("Initial reuse vec dimension mismatch");
    }
        
     LDPC_DE_LUT  de( ens,
                      Nq_Cha,
                      Nq_Msg_vec,
                      maxiter, 
                      var_luts,
                      chk_luts,
                      reuse_vec ,
                      THR_PREC,
                      PE_MAX,
                      LDPC_DE_LUT::ARI,
                      MAX_BISEC_ITER,
                      MAX_LLR_MAGNITUDE,
                      NQ_FINE,
                      lut_table_design);


    int jj = 0;
    vec Pe_min_vec(maxiter);
    ivec iter_vec(maxiter);
    int init_reuse = sum(to_ivec(reuse_vec));
    int num_reuse = maxiter- init_reuse - max_reuse_stages;
    cout << "Starting optimization. Initial reuse stages = " << init_reuse 
         << ", target number of stages = " << max_reuse_stages 
         << ", stages being added = " << num_reuse << endl;

     while(jj<num_reuse){
        // Reset
        Pe_min_vec = ones(maxiter);
        iter_vec = maxiter*ones_i(maxiter);
        // Create threads
        std::vector<std::thread> threads; 
        for(int ii=1; ii<maxiter; ii++){
            if(reuse_vec(ii)==1)    continue;
            else{
              bvec reuse_vec_tmp = reuse_vec;
              reuse_vec_tmp(ii) = 1;
              de.set_reuse_vec(reuse_vec_tmp);
              threads.push_back(std::thread(evolve_thread, de, thr, Pe_max, &Pe_min_vec(ii), &iter_vec(ii)  ));
            }
        }
        // Wait for threads to finish
        for (auto& th : threads) th.join();

        // Find the number of iterations, after which Pe_max was reached
        
        int min_iters_idx = min_index(iter_vec);
        int min_Pe_idx = min_index(Pe_min_vec);

        if(iter_vec(min_iters_idx)==maxiter){ // Pe requirement not fulfilled
            cout << "Could not reach Pe target, scaling down to thr = " << thr*scale_down << endl;
            thr *= scale_down;
        }
        else{
            reuse_vec(min_Pe_idx) = 1;
            cout << "Reached Pe target within " << iter_vec(min_Pe_idx) << " iterations." << endl; 
            jj++;
            cout << "Reuse stage " << jj << ", Adding idx = " << min_Pe_idx << endl
                 << "reuse_vec = " << reuse_vec << endl;
        }
    }
    
    cout << "Finished." << endl << "reuse_vec = " << reuse_vec << endl;
    return EXIT_SUCCESS;
}
