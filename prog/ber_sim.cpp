/*!
 * \file
 * \brief Program to perform LDPC Bit Error Rate Monte Carlo Simulations over an BIAWGN Channel
 * \author Michael Meidlinger
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2016 Michael Meidlinger - All Rights Reserved
 *
 */

#include "LDPC_DE.hpp"
#include "LDPC_Code_LUT.hpp"
#include "LDPC_BER_Sim.hpp"
#include <boost/program_options.hpp>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/optional.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

namespace po = boost::program_options;

using namespace itpp;
using namespace std;
int main(int argc, char **argv){
    
    //=== Parse input parameters
    int seed;
    std::string base_dir;
    std::string custom_name_ext;
    boost::filesystem::path base_dir_path;
    std::string params_filename;
    boost::filesystem::path params_filename_path;
    
    try {
        
        po::options_description desc(
                    "OPTIONS:");
        desc.add_options()
        ("basedir,b", po::value<std::string>(&base_dir)->default_value(boost::filesystem::current_path().string() ), 
	 		"paths in params files are relative to this directory. Default: current direcroy" )
	("custom-name,c", po::value<std::string>(&custom_name_ext)->default_value(""), 
	 		"append this string at the end of the results file name" )
        ("help,h", "produce help message")
        ("params,p", po::value<std::string>(&params_filename), "input parameter file")
        ("seed,s", po::value<int>(&seed)->default_value(0), "random seed")
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
        
        
        base_dir_path  = base_dir;
        if(vm.count("basedir")>1){
            cout << "The base directory can only be specified once" << "\n";
            return 0;
        }
        if (base_dir_path.is_relative()) {
            cout << "Base directory must be specified as absolut path" << "\n";
            return 0;
        }
	
	// Parse custom file name options
	if(vm.count("custom-name")>1){
            cout << "Custom Name can only be specified once" << "\n";
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

        
    //=== Load Simulation
    LDPC_BER_Sim* sim;
    boost::property_tree::ptree param_tree;

    
    if(params_filename_path.is_relative()){
        params_filename_path = base_dir_path / params_filename_path;
    }

    if(!exists(params_filename_path)){
        it_error("Parameter file" << params_filename_path.string() << " does not exist!");
    }
    boost::property_tree::ini_parser::read_ini(params_filename_path.string(), param_tree);
    
    // use this ptree node the check wether the corresponding property nodes exist
    boost::optional<boost::property_tree::ptree& > child;
    
    // If a codec is to be loaded rather than generated, we need to know of which type it is
    std::string codec_type = param_tree.get("Sim.codec_type", "none");
    
    if( (child=param_tree.get_child_optional( "LUT" )) || codec_type == "LUT"  ){
        sim = new LDPC_BER_Sim_LUT(params_filename_path, base_dir_path);
    }
    else if( (child = param_tree.get_child_optional( "BP" )) || codec_type == "BP"){
        sim =  new LDPC_BER_Sim(params_filename_path, base_dir_path);
    }
    else{
        it_error("You must specify the type of decoder in the params file by either \n\
                    * explicitly setting at least one parameter in the [LUT] or [BP] section of the parameter file\
                    * loading a codec via its filename and specifying its type");
	sim = new LDPC_BER_Sim(); //will never be executed, surpress compiler warnings
    }

    //=== Run Simulation and save files
    sim->rand_seed = seed;
    sim->append_custom_name(custom_name_ext);
    sim->load();
    sim->run();
    sim->save();
    
    delete sim;

    
    return EXIT_SUCCESS;
}
