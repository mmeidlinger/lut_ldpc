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

#include "LDPC_Degree_Opt.hpp"
#include <boost/program_options.hpp>

#define MAX_LLR_MAGNITUDE 25
#define MAX_BISEC_ITER 50
#define NQ_FINE 5000
#define MAX_ITERS_DEGREE_OPT 2000

namespace po = boost::program_options;
using namespace std;
using namespace itpp;



int main(int argc, char **argv){

    cout << "Called program via" << endl << argv[0];
    for(int ii=1; ii<argc; ii++){
        cout << " " << argv[ii];
    } 
    cout << endl;

    bool lut_de=false; // true if LUT decoding is used for DE, false if BP decoding is used
    bool min_approx=false;   // true if min-LUT algorithm should be used
    double R_target;    
    int Nq;
    int iters_de;
    string ensemble_filename;
    double scale_down;
    LDPC_Ensemble ens;
    std::vector<string> degree_dist;
    double delta;
    double Pe_max;
    double threshold_precision;
    string lut_tree_design, lut_table_design;
    bool var_opt = false;
    bool chk_opt = false;
    bool lam2constraint = false;

    try {
        po::options_description desc(
                    "OPTIONS:");
        desc.add_options()
        ("lut,L", 
            po::bool_switch(&lut_de),
	 		"Use LUT decoding algorithm" )
        ("min-approx,m",
            po::bool_switch(&min_approx),
	 		"Approximate Check Node Updates" )
        ("quant-bits,q",
            po::value<int>(&Nq),
            "Number of quantization bits")
        ("ensemble,e",
            po::value<string>(&ensemble_filename),
            "Filename for initial ensemble")
        ("iterations,i",
            po::value<int>(&iters_de)->default_value(2000),
            "Number of Message passing iterations")
        ("degree-dist,d",
            po::value<std::vector<string> >(&degree_dist)->multitoken(),
            "Degree ditribution is the form,  \"VN_degrees / VN_probabilities / CN_degrees / CN_probabilities \"")
        ("rate,R",
            po::value<double>(&R_target)->required(),
            "Target rate")
        ("scale-down,s",
            po::value<double>(&scale_down)->default_value(0.995),
            "Scale down threshold by this value if an updated ensemble fails to converge")
        ("pmax,p",
            po::value<double>(&Pe_max)->default_value(1e-16),
            "Convergence error probability")
        ("delta",
            po::value<double>(&delta)->default_value(1e-2),
            "Distance constraint for Linear Program formulation to be valid")
        ("threshold-precision,t",
            po::value<double>(&threshold_precision)->default_value(1e-5),
            "Bisection search precision")
        ("lut-table-design",
            po::value<string>(&lut_table_design)->default_value("joint_root"),
            "Strategy for LUT table creation")
        ("lut-tree-design",
            po::value<string>(&lut_tree_design)->default_value("auto_bin_balanced"),
            "Strategy for LUT table creation")
        ("var-opt,v", 
            po::bool_switch(&var_opt),
            "Optimize variable node degree distribution")
        ("chk-opt,c", 
            po::bool_switch(&chk_opt),
            "Optimize check node degree distribution")
        ("lam2constraint",
            po::bool_switch(&lam2constraint),
            "Constraint maximum value of variable node edge degree 2 according to theoretical limits")
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

        //==== Check optimization metho
        it_assert(var_opt || chk_opt, "Specify what degree distributions to optimize!");
        //==== Parse initial degree distribution 
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
        // Echo initial ensemble
        cout << "Successfully read ensemble with rate " <<  ens.get_rate() << endl << ens << endl;
        // Check rate
        double R_init = ens.get_rate();
        it_assert(R_init < R_target, "The rate of the initial ensemble (=" << R_init << ") must be smaller than the rate target!" );
        
        //==== Parse resolution
        if (vm.count("quant-bits")==0){
            if(lut_de)  Nq = 4;
            else        Nq = 11;
        }
        
        //==== check for exclusive options 
        if(lut_de==false && vm["lut-tree-design"].defaulted()==false ){
            it_warning("--lut-tree-design only applies if -L has been set, ignoring.");
        }
        if(!lut_de==false && vm["lut-table-design"].defaulted()==false ){
            it_warning("--lut-table-design only applies if -L has been set, ignoring.");
        }


    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

   


    LDPC_DE* de; 
    if(lut_de == false){
        if(min_approx) {it_error("Min sum density evolution not implemented yet!"); }
        de = new LDPC_DE_BP(ens,  Nq, MAX_LLR_MAGNITUDE);
    }
    else{
        int Nq_Cha = pow2i(Nq);
        ivec Nq_Msg_vec =   pow2i(Nq)*ones_i(iters_de);
        bvec reuse_vec = zeros_b(iters_de);
        Array<Array<LUT_Tree>> var_luts;
        Array<Array<LUT_Tree>> chk_luts;
        
        get_lut_tree_templates(lut_tree_design, ens, Nq_Msg_vec, Nq_Cha, min_approx, var_luts, chk_luts);
        
        
        de = new LDPC_DE_LUT( ens,
                              Nq_Cha,
                              Nq_Msg_vec,
                              iters_de,
                              var_luts,
                              chk_luts,
                              reuse_vec ,
                              threshold_precision ,
                              Pe_max,
                              LDPC_DE_LUT::ARI,
                              MAX_BISEC_ITER,
                              MAX_LLR_MAGNITUDE,
                              NQ_FINE,
                              lut_table_design);
    }
    
    LDPC_Degree_Opt_LP deg_opt(de,
                               R_target,
                               var_opt,
                               chk_opt,
                               delta,
                               scale_down, 
                               lam2constraint,
                               MAX_ITERS_DEGREE_OPT);
    
    deg_opt.run(ens);
    
    
    return EXIT_SUCCESS;
}
