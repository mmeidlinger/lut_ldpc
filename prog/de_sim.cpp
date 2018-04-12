/*!
 * \file
 * \brief Program to perform multithreaded Density Evolution simulations
 * \author Michael Meidlinger
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2017 Michael Meidlinger - All Rights Reserved
 *
 * This file is part of lut_ldpc, a software suite for simulating and designing
 * LDPC decodes based on discrete Lookup Table (LUT) message passing
 *
 * lut_ldpc is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * lut_ldpc distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with lut_ldpc.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include "LDPC_DE.hpp"

#include <boost/program_options.hpp>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/optional.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <thread>



using namespace itpp;
using namespace std;
using namespace lut_ldpc;

// Declarations
void de_lut_thread_call(LDPC_DE_LUT de, double* thr, int* iters);
void de_bp_thread_call(LDPC_DE_BP de, double* thr, int* iters);

//! Perform Density evolution for a LUT decoder according to parameters
void de_sim_lut(const boost::property_tree::ptree& param_tree);
//! Perform Density evolution for a BP decoder according to parameters
void de_sim_bp(const boost::property_tree::ptree& param_tree);


namespace po = boost::program_options;

int main(int argc, char **argv)
{
    //=== Parse input parameters
    std::string params_filename;
    boost::filesystem::path params_filename_path;
    
    try {
        
        po::options_description desc(
                                     "OPTIONS:");
        desc.add_options()
        ("help,h", "produce help message")
        ("params,p", po::value<std::string>(&params_filename), "input parameter file")
        ;
        
        
        
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        
        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }
        
        params_filename_path = params_filename;
        if (vm.count("params")==0) {
            cout << "No input parameters specified. To learn more, use the --help option." << "\n";
            return 0;
        }
        
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
    
    //=== Parse parameter file
    boost::property_tree::ptree param_tree;
    boost::property_tree::ini_parser::read_ini(params_filename, param_tree);
    boost::optional<std::string> input_string;
    // use this ptree node the check wether the corresponding property nodes exist
    boost::optional<boost::property_tree::ptree&> child;
    
    
    if( (child=param_tree.get_child_optional( "LUT" )) ){
        de_sim_lut(param_tree);
    }
    else if( (child = param_tree.get_child_optional( "BP" ))){
        de_sim_bp(param_tree);
    }
    else{
        it_error("You must specify the type of decoder in the params file by explicitly setting at least one parameter in \
                 either the [LUT] or [BP] section of the parameter file!");
    }
    
    
    

    return EXIT_SUCCESS;
}


void de_lut_thread_call(LDPC_DE_LUT de, double* thr, int* iters){
    std::cout << "Starting thread" << std::endl;
    *iters = de.bisec_search(*thr);
    std::cout << "Done with thread" << std::endl;
}

void de_bp_thread_call(LDPC_DE_BP de, double* thr, int* iters){
    std::cout << "Starting thread" << std::endl;
    *iters = de.bisec_search(*thr);
    std::cout << "Done with thread" << std::endl;
}

void de_sim_lut(const boost::property_tree::ptree& param_tree)
{
    // Load the ensemble an print out rate
    std::string ensemble_filename = param_tree.get<std::string>("Sim.ensemble_filename");
    LDPC_Ensemble ens(ensemble_filename);
    
    cout << "Density evolution simualtion for ensemble of Rate " << ens.get_rate() << endl;
    
    // Simulation settings
    double thr_min   = param_tree.get("Sim.thr_min", 1e-9);
    boost::optional<double> thr_max   =  param_tree.get_optional<double>("Sim.thr_max");
    if(!thr_max)    thr_max = rate_to_shannon_thr(ens.get_rate());
    double thr_prec   = param_tree.get("Sim.thr_prec", 1e-4);
    double Pe_max    =  param_tree.get("Sim.Pe_max", 1e-9);
    ivec maxiter_de_vec     =  param_tree.get("Sim.maxiter_de", "1000");
    int maxiter_bisec     =  param_tree.get("Sim.maxiter_bisec", 50);
    int max_ni_de_iters =  param_tree.get("Sim.max_ni_de_iters", 30);
    
    double LLR_max    =  param_tree.get("Sim.LLR_max", 25.0);
    std::string results_name = param_tree.get<std::string>("Sim.results_name");
    
    //LUT-specific settings
    imat qbits = param_tree.get("LUT.qbits",  "3 3; 4 4");
    ivec Nq_msg_vec = param_tree.get("LUT.Nq_msg_vec",  " ");
    ivec reuse_iter_vec     =  param_tree.get("LUT.reuse_iter_vec", "0");    
    bvec reuse_vec_in =  param_tree.get("LUT.reuse_vec", "");      
    bool  min_lut = param_tree.get("LUT.min_lut", false);
    string tree_mode = param_tree.get("LUT.tree_mode", "auto_bin_balanced");
    int  Nq_fine = param_tree.get("LUT.Nq_fine", 5000);
    std::string irregular_design_strategy = param_tree.get("LUT.irregular_design_strategy", "joint_root");
    
    
    // Start the evolution
    int num_threads;
    if( reuse_iter_vec.size() == 1 && qbits.rows()==1 && maxiter_de_vec.size()>=1 ){
        num_threads = maxiter_de_vec.size();
    }
    else if( reuse_iter_vec.size() == 1  && maxiter_de_vec.size()==1 && qbits.rows()>=1){
        num_threads = qbits.rows();
    }
    else if( reuse_iter_vec.size() >= 1  && maxiter_de_vec.size()==1 && qbits.rows()==1){
        num_threads = reuse_iter_vec.size();
    }
    else{
        it_error("de_sim_lut(): This function either sweeps over LUT resolutions, number of iterations or reuse factors!");
	num_threads = 0; // Won't be executed, surpress compiler warnings
    }
    
    double* thresholds = new double [num_threads];
    int* bisec_iters = new int [num_threads];
    std::vector<std::thread> threads;
    Array<LDPC_DE_LUT> de_array(num_threads);
    
    
    for(int nn=0; nn< num_threads; nn++){
        // Get trees
        Array<Array<LUT_Tree> > var_luts;
        Array<Array<LUT_Tree> > chk_luts;
        
        int Nq_msg=0;
        int Nq_cha=0;
        int maxiter_de=0;
        int reuse_iters=0;
        
        
        if(maxiter_de_vec.size() == num_threads){
            Nq_msg = qbits(0,1);
            Nq_cha = qbits(0,0);
            reuse_iters = reuse_iter_vec(0);
            maxiter_de = maxiter_de_vec(nn);
        }
        else if(qbits.rows() == num_threads){
            Nq_msg = qbits(nn,1);
            Nq_cha = qbits(nn,0);
            reuse_iters = reuse_iter_vec(0);
            maxiter_de = maxiter_de_vec(0);
        }
        else if(reuse_iter_vec.size() == num_threads){
            Nq_msg = qbits(0,1);
            Nq_cha = qbits(0,0);
            reuse_iters = reuse_iter_vec(nn);
            maxiter_de = maxiter_de_vec(0);
        }
        else{
            it_error("de_sim_lut(): This function either sweeps over LUT resolutions, number of iterations pre reuse factors!");
	        maxiter_de = 0; // Won't be executed, surpress compiler warnings
            reuse_iters = 0; // Won't be executed, surpress compiler warnings
        }
       
        get_lut_tree_templates(tree_mode, ens, ones_i(maxiter_de)*pow2i(Nq_msg), pow2i(Nq_cha), min_lut, var_luts, chk_luts);
        
        // Build reuse vector
        bvec reuse_vec;
        if(reuse_vec_in.length()==0){
            reuse_vec = zeros_b(maxiter_de);
            it_assert(reuse_iters>=0, "de_sim_lut(): thread " << nn <<": Reuse iters must be >= 0");
            int reuse_iters_tmp = 0;
            for(int ii=1; ii< maxiter_de-1; ii++){
                if(reuse_iters_tmp < reuse_iters){
                    reuse_vec(ii)=true;
                    reuse_iters_tmp++;
                }
                else{
                    reuse_vec(ii)=false;
                    reuse_iters_tmp=0;
                }
            } 
        }
        else{
           reuse_vec = reuse_vec_in; 
        }
        //reuse_vec = zeros_b(maxiter_de);
        //reuse_vec((int)maxiter_de/2)=true;

        //cout << "Reuse vec = " << reuse_vec << endl;
        
        ivec Nq_msg_vec_tmp;
        if(Nq_msg_vec.size() == maxiter_de) // If a message vec was explicitly set, we use that
            Nq_msg_vec_tmp = to_ivec(pow2(Nq_msg_vec));
        else
            Nq_msg_vec_tmp = ones_i(maxiter_de)*pow2i(Nq_msg);
        LDPC_DE_LUT de(   ens,
                          pow2i(Nq_cha),
                          Nq_msg_vec_tmp,
                          maxiter_de,
                          var_luts,
                          chk_luts,
                          reuse_vec,
                          thr_prec ,
                          Pe_max,
                          LDPC_DE::ARI,
                          maxiter_bisec,
                          LLR_max,
                          Nq_fine,
                          irregular_design_strategy);
            
            de.set_bisec_window(thr_min, *thr_max);
            de.set_exit_conditions(maxiter_de, maxiter_bisec, max_ni_de_iters, Pe_max, thr_prec);
            threads.push_back(std::thread(de_lut_thread_call, de, &thresholds[nn], &bisec_iters[nn]  ));
            de_array(nn)=de;

    }
    
    // Wait for threads to finish
    for (auto& th : threads) th.join();
    
    // Write results to vectors and delete threads
    ivec bisec_iters_vec(num_threads);
    vec threshold_vec(num_threads);
    for(int nn=0; nn< num_threads; nn++){
        bisec_iters_vec(nn) = bisec_iters[nn];
        threshold_vec(nn) = thresholds[nn];
    }
    // clear memory
    delete [] thresholds;
    delete [] bisec_iters;
    
    // Calculate stable degrees for thresholds
    vec lam2stable_vec = zeros(threshold_vec.length());
    vec rho = ens.get_chk_degree_dist();
    for(int nn=0; nn< num_threads; nn++){
        int Nq_msg=0;
        int Nq_cha=0;
        if(maxiter_de_vec.size() == num_threads){
            Nq_msg = qbits(0,1);
            Nq_cha = qbits(0,0);
        }
        else if(qbits.rows() == num_threads){
            Nq_msg = qbits(nn,1);
            Nq_cha = qbits(nn,0);
        }
        else if(reuse_iter_vec.size() == num_threads){
            Nq_msg = qbits(0,1);
            Nq_cha = qbits(0,0);
        }
        else{
            it_error("de_sim_lut(): This function either sweeps over LUT resolutions, number of iterations pre reuse factors!");
        }
        lam2stable_vec(nn) = get_lam2stable_lut(threshold_vec(nn), rho ,  pow2(Nq_cha), pow2(Nq_msg) , 25.0, 5000);
    }
     
    std::ofstream file;
    
    file.open(results_name.c_str());
    it_assert(file.is_open(),
              "de_sim_bp(): Could not open file \""
              << results_name << "\" for writing");
    file << "==== DE Threshold for ensemble file " << ensemble_filename
    << " (Rate = " << ens.get_rate() << ", BI-AWGN channel) " << endl
    << "  Active Variable node degrees: " << ens.sget_degree_lam() << endl
    << "  pmf of Variable node edges: " << ens.sget_lam() << endl
    << "  Active Check node degrees: " << ens.sget_degree_rho() << endl
    << "  pmf of Check node edges: " << ens.sget_rho() << endl;
    
    file << "-- SIMULATION PARAMETERS"
    << "  Search Window = [" << thr_min << ", " << *thr_max << "]" << endl
    << "  Threshold precision = " << thr_prec << endl
    << "  Convergence error probability = " << Pe_max << endl
    << "  Maximum Number of message passing iterations = " << maxiter_de_vec << endl
    << "  MinLut Algorithm used = " << min_lut << endl
    << "  LUT Tree design mode = " << tree_mode << endl
    << "  LUT table design mode = " << irregular_design_strategy << endl
    << "  LUT reuse iter vec = " << reuse_iter_vec << endl
    << "  Non improving iterations tolerated before terminating = " << max_ni_de_iters << endl
    << "  Resolutions [channel bits, message bits; ...] = " << qbits <<  endl
    << "  Program git version = " << gitversion << endl
    << "  Bisection iterations until convergence = "  << bisec_iters_vec << endl
    << "  Stable lam2 degrees at thresholds = " << lam2stable_vec << endl
    << "  Threshold(s) found = " << threshold_vec << endl
    << "  Eb/N0 corresponding to thresholds = " << sig2snr(ens.get_rate(), threshold_vec) << endl;
    
    
    if(Nq_msg_vec.size() > 0){
        file << "  Unique Nq_msg_vec = " << Nq_msg_vec << endl;
    }
    file << endl;
    

    // Get Pe trace at threshold
    vec Pe_trace;
    if(num_threads == 1){
        cout << "Calculating Pe trace for threshold " << threshold_vec(0) << endl;
        mat dummy_mat;        
        de_array(0).evolve(threshold_vec(0), true, false, dummy_mat, Pe_trace);
        // The following two lines are for testing without having an explicit program to interact
        // bvec reuse_vec = de_array(0).evolve_adaptive_reuse(threshold_vec(0)*0.99, Pe_trace, 0.1, 1e-5, 13);
        // cout << reuse_vec << endl;
        file << "  Pe_trace = "<< Pe_trace << endl;
    }

    file.close();

}

void de_sim_bp(const boost::property_tree::ptree& param_tree)
{
    // Load the ensemble
    std::string ensemble_filename = param_tree.get<std::string>("Sim.ensemble_filename");
    LDPC_Ensemble ens(ensemble_filename);
    
    cout << "Density evolution simualtion for ensemble of Rate " << ens.get_rate() << endl;
    
    // Simulation settings
    double thr_min   = param_tree.get("Sim.thr_min", 1e-9);
    boost::optional<double> thr_max   =  param_tree.get_optional<double>("Sim.thr_max");
    if(!thr_max)    thr_max = rate_to_shannon_thr(ens.get_rate());
    double thr_prec   = param_tree.get("Sim.thr_prec", 1e-4);
    double Pe_max    =  param_tree.get("Sim.Pe_max", 1e-9);
    ivec maxiter_de_vec     =  param_tree.get("Sim.maxiter_de", "1000");
    int maxiter_bisec     =  param_tree.get("Sim.maxiter_bisec", 50);
    int max_ni_de_iters =  param_tree.get("Sim.max_ni_de_iters", 5);
    double LLR_max    =  param_tree.get("Sim.LLR_max", 25.0);
    std::string results_name = param_tree.get<std::string>("Sim.results_name");
    
    
    //BP-specific settings
    int Nq = param_tree.get("BP.qbits", 10);
    bool  min_sum = param_tree.get("BP.min_sum", false);
    it_assert(min_sum == false, "de_sim_bp(): Min sum density evolution not implemented yet!");
    
    
    boost::optional<string> tree_file   =  param_tree.get_optional<string>("LUT.trees_filename");
    if(!tree_file){
        string tree_mode = param_tree.get("LUT.tree_mode", "auto_bin_balanced");
    }
    else{
        
    }
    
    int num_threads = maxiter_de_vec.size();
    double* thresholds = new double [num_threads];
    int* bisec_iters = new int [num_threads];
    std::vector<std::thread> threads;

    
    for(int nn=0; nn< num_threads; nn++){
        // Start the evolution
        LDPC_DE_BP de(ens, Nq, LLR_max);
        de.set_bisec_window(thr_min, *thr_max);
        de.set_exit_conditions(maxiter_de_vec(nn), maxiter_bisec, max_ni_de_iters, Pe_max, thr_prec);
        threads.push_back(std::thread(de_bp_thread_call, de, &thresholds[nn], &bisec_iters[nn]  ));
        
    }
    
    // Wait for threads to finish
    for (auto& th : threads) th.join();
    
    // Write results to vectors and delete threads
    ivec bisec_iters_vec(num_threads);
    vec threshold_vec(num_threads);
    for(int nn=0; nn< num_threads; nn++){
        bisec_iters_vec(nn) = bisec_iters[nn];
        threshold_vec(nn) = thresholds[nn];
    }
    // clear memory
    delete [] thresholds;
    delete [] bisec_iters;
    
    // Calculate stable degrees for thresholds
    vec lam2stable_vec = zeros(threshold_vec.length());
    vec rho = ens.get_chk_degree_dist();
    for(int nn=0; nn< num_threads; nn++){
        lam2stable_vec(nn) = get_lam2stable_cbp(threshold_vec(nn), rho);
    }
    
    std::ofstream file;

    file.open(results_name.c_str());
    it_assert(file.is_open(),
              "de_sim_bp(): Could not open file \""
              << results_name << "\" for writing");
    file << "==== DE Threshold for ensemble file " << ensemble_filename
         << " (Rate = " << ens.get_rate() << ", BI-AWGN channel) " << endl
         << "  Active Variable node degrees: " << ens.sget_degree_lam() << endl
         << "  pmf of Variable node edges: " << ens.sget_lam() << endl
         << "  Active Check node degrees: " << ens.sget_degree_rho() << endl
         << "  pmf of Check node edges: " << ens.sget_rho() << endl;
    
    file << "-- SIMULATION PARAMETERS" << endl
         << "  Search Window = [" << thr_min << ", " << *thr_max << "]" << endl
         << "  Threshold precision = " << thr_prec << endl
         << "  Convergence error probability = " << Pe_max << endl
         << "  Maximum Number of message passing iterations = " << maxiter_de_vec << endl
         << "  MinSum Approximation used = " << min_sum << endl
         << "  Non improving iterations tolerated before terminating = " << max_ni_de_iters << endl
         << "  Resolution of discrete pmfs = " << Nq << " bit" << endl
         << "  Maximum LLR magnitude = " << LLR_max << endl
         << "  Program git version = " << gitversion << endl
         << "  Bisection iterations until convergence = " << bisec_iters_vec << endl
         << "  Stable lam2 degrees at thresholds = " << lam2stable_vec << endl
         << "  Threshold(s) found = " << threshold_vec << endl
         << "  Eb/N0 corresponding to thresholds = " << sig2snr(ens.get_rate(), threshold_vec) << endl << endl;
    file.close();
    
}
