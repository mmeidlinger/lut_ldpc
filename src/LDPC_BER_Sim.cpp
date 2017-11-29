/*!
 * \file LDPC_BER_Sim.cpp
 * \brief Implementation of Bit Error Rate simulations (BER) for LDPC codes over and Binary-Input AWGN channel
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

#include "LDPC_BER_Sim.hpp"

using namespace lut_ldpc;
using namespace std;
using namespace itpp;
using namespace boost::filesystem;



// ----------------------------------------------------------------------
// LDPC_BER_Sim
// ----------------------------------------------------------------------

LDPC_BER_Sim::LDPC_BER_Sim(const boost::filesystem::path& params_file_path, const boost::filesystem::path& base_dir_path)
{
    boost::property_tree::ptree param_tree;
    if(!exists(params_file_path)){
        it_error("Parameter file" << params_file_path.string() << " does not exist!");
    }
    this->params_file_path = params_file_path;
    
    boost::property_tree::ini_parser::read_ini(params_file_path.string(), param_tree);
    
    // Simulation settings
    SNRdB   = param_tree.get<std::string>("Sim.SNRdB");
    Nframes = param_tree.get("Sim.Nframes", 1e5);
    Nfers   = param_tree.get("Sim.Nfers", 20);
    ber_min   = param_tree.get("Sim.ber_min", 1e-7);
    fer_min   = param_tree.get("Sim.fer_min", 1e-5);
    rand_seed_offset = param_tree.get("Sim.rand_seed_offset", 0);
    save_codec = param_tree.get("Sim.save_codec", 0);
    custom_name = param_tree.get("Sim.custom_name", "");
    results_prefix = param_tree.get("Sim.results_prefix", "RES");
    results_dir = param_tree.get("Sim.results_dir", "results");
    codes_dir = param_tree.get("Sim.codes_dir", "codes");
    codec_filename = param_tree.get("Sim.codec_filename", "");
    
    
    
    // LDPC Code Settings
    parity_filename = param_tree.get("LDPC.parity_filename", "");
    zero_codeword = param_tree.get("LDPC.zero_codeword", true);
    save_permuted = param_tree.get("LDPC.save_permuted", false);
    parity_check_iter = param_tree.get("LDPC.parity_check_iter", true);
    
    // Settings for Belief Propagation decoder
    max_iter =param_tree.get("BP.max_iter", 30);
    llr_calc_d1 = param_tree.get("BP.qllr_scale_res" , 12);
    llr_calc_d2 = param_tree.get("BP.qllr_table_size" , 300);
    llr_calc_d3 = param_tree.get("BP.qllr_spacing_res" , 7);
    llr_calc_d4 = param_tree.get("BP.qllr_total_res" , 8*(int)sizeof(QLLR)-4);
    
    
    // Set path objects and create folders if necessarry
    boost::system::error_code mkdir_error;
    
    codes_path =  codes_dir;
    if(codes_path.is_relative()){
        codes_path = base_dir_path / codes_path;
    }
    create_directories( codes_path, mkdir_error );
    if(mkdir_error){
        it_error("LDPC_BER_Sim::LDPC_BER_Sim(): Could not create directory" << codes_path.string());
    }
    
    results_path = results_dir;
    if(results_path.is_relative()){
        results_path = base_dir_path / results_path;
    }
    create_directories( results_path, mkdir_error );
    if(mkdir_error){
        it_error("LDPC_BER_Sim::LDPC_BER_Sim(): Could not create directory" << results_path.string());
    }
    
}

std::string LDPC_BER_Sim::gen_filename() const {
    it_assert(decoder_set, "LDPC_BER_Sim::gen_filename(): Decoder has not been set!");
    std::stringstream fn;
    fn  << results_prefix
        << "_N" << codeword_length
        << "_R" << C->get_rate()
        << "_maxIter" << max_iter
        << "_zcw" << (int) zero_codeword
        << "_frames" << Nframes
        << custom_name;
    return fn.str();
}

void LDPC_BER_Sim::append_custom_name(std::string& custom_name_ext){
    custom_name += custom_name_ext;
}

void LDPC_BER_Sim::run(){
    it_assert(decoder_set, "LDPC_BER_Sim::run(): Decoder has not been set!");

    // Initialize clock to save run time
    std::clock_t start;
    double duration;
    start = std::clock();

    // Initialize random seed
    RNG_reset(rand_seed+rand_seed_offset);
    
    // Run Simulations
    int ss = 0;
    while(ss< length(SNRdB))
    {
        bool exit_cond = sim_snr_point(SNRdB(ss));
        ss++;
        if(exit_cond == true)
            break;
    }
    
    // Add additional empty SNR points if neccessarry
    while(ss< length(SNRdB)){
        BERC berc_data;
        BERC berc_uncoded;
        BLERC ferc;
        ferc.set_blocksize(dataword_length);
        results.add_snr_point(SNRdB(ss), berc_data, berc_uncoded, ferc);
        ss++;
    }
    
    // Stop run time
    double runtime = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    results.save_runtime(runtime);
    cout << "Done simulating. Runtime = " << runtime << " seconds" << endl;
}

void LDPC_BER_Sim::load() {

    if(codec_filename.length() == 0){
        path parity_path = codes_path / (parity_filename + ".alist");
        if(!exists(parity_path)){
            it_error("Parity file" << parity_path.string() << " does not exist!");
        }
        LDPC_Parity* H = new LDPC_Parity(parity_path.string());
       
        LDPC_Code* C_ldpc;
        // Load Generator if needed and present
        if(zero_codeword == false){
            LDPC_Generator_Systematic* G;
            path generator_path = codes_path / (parity_filename + ".gen.it");
            if(exists(generator_path)){
                G = new LDPC_Generator_Systematic();
                G->load(generator_path.string());
            }
            else{
                it_info_debug("No generator file found" << std::endl
                              << "Calculating Generator based on Parity Check Matrix ..." << std::endl);
                G = new LDPC_Generator_Systematic(H);
                
                // Save newly generated Generator and correspondingly permuted Parity
                if(save_permuted){
                    H->save_alist(parity_path.string());
                    it_file f(generator_path.string());
                    f<< Name("Fileversion") << 2;
                    f.close();
                    G->save(generator_path.string());
                    it_info_debug("done");
                }
            }
            // Create Codec capable of encoding and decoding
            C_ldpc = new LDPC_Code(H, G);
            decoder_set = true;
            encoder_set = true;
        }
        else{
            // Create Codec capable of decoding only
            C_ldpc = new LDPC_Code(H);
            decoder_set = true;
        }
        C_ldpc->set_exit_conditions(max_iter,parity_check_iter,parity_check_iter);
        C_ldpc->set_llrcalc(LLR_calc_unit(llr_calc_d1, llr_calc_d2, llr_calc_d3, llr_calc_d4));
        
        codeword_length = C_ldpc->get_nvar();
        dataword_length = C_ldpc->get_ninfo();
        C = C_ldpc;
        
        // Save codec
        if(rand_seed == save_codec){
            path results_path_subdir = results_path / gen_filename();
            boost::system::error_code mkdir_error;
            create_directories( results_path_subdir, mkdir_error );
            if(mkdir_error){
                it_error("LDPC_BER_Sim::load(): Could not create directory" << results_path_subdir.string());
            }
            C_ldpc->save_code((results_path_subdir / "bp_codec.it").string());
        }
    }
    else{
        path codec_path = codes_path / (codec_filename);
        if(!exists(codec_path)){
            it_error("Codec file" << codec_path.string() << " does not exist!");
        }
        LDPC_Generator_Systematic* G  = new LDPC_Generator_Systematic();
        LDPC_Code* C_ldpc = new LDPC_Code(codec_path.string(), G);
        C = C_ldpc;
        C_ldpc->set_exit_conditions(max_iter,parity_check_iter,parity_check_iter);
        C_ldpc->set_llrcalc(LLR_calc_unit(llr_calc_d1, llr_calc_d2, llr_calc_d3, llr_calc_d4));
        decoder_set = true;
        if (G->is_initialized() == true){
            encoder_set = true;
        }
        else{
            delete G;
            G = nullptr;
        }
        codeword_length = C_ldpc->get_nvar();
        dataword_length = C_ldpc->get_ninfo();
        max_iter = C_ldpc->get_nrof_iterations();

    }
    
    // Create results
    results = LDPC_BER_Sim_Results(codeword_length, codeword_length - dataword_length);
}

bool LDPC_BER_Sim::sim_snr_point(double snr) {
    // Noise variance is N0/2 per dimension
    double N0 = itpp::pow10( -snr / 10.0) / C->get_rate();
    AWGN_Channel chan(N0 / 2);
    BERC berc_data;         // Counters for BER of data bits
    BERC berc_uncoded;      // Counters for BER if hard output demodulation would be used
    BLERC ferc;         // Counter for coded FER
    
    BPSK Mod;
    
    int N = codeword_length;
    int K = dataword_length;
    
    ferc.set_blocksize(K);
    for (int64_t ff = 0; ff < Nframes; ff++) {
        bvec bitsin, bitsout, bitsin_coded, bitsout_slicer ;
        // Generate Transmit Signal
        if(zero_codeword){
            bitsin = zeros_b(K);
            bitsin_coded = zeros_b(N);
        }
        else if(encoder_set){
            bitsin = randb(K);
            bitsin_coded = C->encode(bitsin);
        }
        else{
            it_error("Non zero codewords require the encoder to be set!");
        }
        vec s = Mod.modulate_bits(bitsin_coded);
        // Received data
        vec x = chan(s);
        // Demodulate
        vec softbits = Mod.demodulate_soft_bits(x, N0);
        bitsout_slicer = Mod.demodulate_bits(x);
        // Decode the received bits
        bitsout = C->decode(softbits);
        
        // Count the number of errors
        berc_data.count(bitsin, bitsout);
        berc_uncoded.count(bitsin_coded,bitsout_slicer);
        ferc.count(bitsin, bitsout);
        
        // Break if maximum number of frame errors is reached
        if (ferc.get_errors() > Nfers) break;
        
    }
    
    // Output
    std::cout << "SNR = " << snr << "  Simulated "
    << ferc.get_total_blocks() << " frames and "
    << berc_data.get_total_bits() << " data bits. "
    << "Obtained " << berc_data.get_errors() << " data bit errors. "
    << " Data BER: " << berc_data.get_errorrate()
    << " Uncoded BER: " << berc_uncoded.get_errorrate()
    << " FER: " << ferc.get_errorrate() << std::endl << std::flush;
    
    
    // Write SNR Point
    results.add_snr_point(snr, berc_data, berc_uncoded, ferc);
    
    // Break if Minimum data BER is reached
    if (berc_data.get_errorrate() < ber_min || ferc.get_errorrate() < fer_min )
        return true;
    else
        return false;
}

// ----------------------------------------------------------------------
// BER_Sim_Results
// ----------------------------------------------------------------------

void LDPC_BER_Sim::save(){
    path results_path_subdir = results_path / gen_filename();
    boost::system::error_code mkdir_error;
    create_directories( results_path_subdir, mkdir_error );
    if(mkdir_error){
        it_error("LDPC_BER_Sim::save(): Could not create directory" << results_path_subdir.string());
    }
    
    // Build full filename
    stringstream fstream;
    fstream << gen_filename() << "_rseed"<< setfill('0') << setw(4) << rand_seed + rand_seed_offset << ".it";
    // Call OS and Path agnostic function to write the results file
    results.write_itfile( (results_path_subdir / fstream.str() ).string());
    // Copy params file into results directory
    try {
        copy_file(params_file_path, results_path_subdir / params_file_path.filename() , copy_option::fail_if_exists);
    }
    catch (const boost::system::system_error &err) {
        // Only throw an error if the error is not due to the file existing already
        if (err.code() != (boost::system::errc::file_exists))
            throw;
    }
    
}

void LDPC_BER_Sim_Results::write_itfile(const std::string& filename){
    
    it_file f;
    f.open(filename);
    // Write Simulation Results
    f << Name("sim_SNRdB") << sim_SNRdB;
    f << Name("sim_Nframes") << to_vec(sim_Nframes);
    f << Name("sim_Ndatabits") << to_vec(sim_Ndatabits);
    f << Name("sim_frame_errors") << to_vec(sim_frame_errors);
    f << Name("sim_data_bit_errors") << to_vec(sim_data_bit_errors);
    f << Name("sim_uncoded_bit_errors") << to_vec(sim_uncoded_bit_errors);
    // Write dependent LDCP params needed to interpret the results in a meaningful way
    f << Name("ldpc_nvar") << to_vec(ldpc_nvar);
    f << Name("ldpc_nchk") << to_vec(ldpc_nchk);
    f << Name("ldpc_code_rate") << to_vec(ldpc_code_rate);
    f << Name("runtime") << runtime;
    f << Name("gitversion") << string(gitversion);
    f.close();
    
    
}
void LDPC_BER_Sim_Results::add_snr_point(double snr, const BERC& berc_data, const BERC& berc_uncoded, const BLERC& ferc){
    sim_SNRdB = concat(sim_SNRdB, snr);
    sim_Nframes = concat(sim_Nframes, (int64_t)ferc.get_total_blocks());
    sim_Ndatabits = concat(sim_Ndatabits, (int64_t)berc_data.get_total_bits());
    sim_frame_errors = concat(sim_frame_errors, (int64_t)ferc.get_errors());
    sim_data_bit_errors = concat(sim_data_bit_errors, (int64_t)berc_data.get_errors());
    sim_uncoded_bit_errors = concat(sim_uncoded_bit_errors, (int64_t)berc_uncoded.get_errors());
}

// ----------------------------------------------------------------------
// LDPC_BER_Sim_LUT
// ----------------------------------------------------------------------

LDPC_BER_Sim_LUT::LDPC_BER_Sim_LUT(const boost::filesystem::path& params_file_path, const boost::filesystem::path& base_dir_path):
    LDPC_BER_Sim(params_file_path, base_dir_path) // call uperclass constructor with same params_file_path
{
    boost::property_tree::ptree param_tree;
    boost::property_tree::ini_parser::read_ini(params_file_path.string(), param_tree);
    boost::optional<std::string> input_string;
    

    int Nqbits_Msg;
    // Settings for LUT decoder
    max_iter = param_tree.get("LUT.max_iter", 30);
    design_thr = param_tree.get_optional<double>("LUT.design_thr");
    design_SNRdB = param_tree.get_optional<double>("LUT.design_SNRdB");
    
    it_assert(design_thr || design_SNRdB, "LDPC_BER_Sim_LUT::LDPC_BER_Sim_LUT(): No design SNR or noise thresold specified");
    
    Nq_Cha   = pow2i( param_tree.get("LUT.qbits_channel", 4) );
    Nqbits_Msg = param_tree.get("LUT.qbits_message_uniform",3);
    
    input_string = param_tree.get_optional<std::string>("LUT.qbits_messages");
    if(input_string)    Nq_Msg = to_ivec(pow2(vec(*input_string)));
    else Nq_Msg = to_ivec(ones(max_iter)*pow2(Nqbits_Msg));
            
    

    
    tree_mode = param_tree.get("LUT.tree_mode", "auto_bin_balanced");
    trees_dir = param_tree.get("LUT.trees_dir", "trees");
    trees_filename = param_tree.get("LUT.trees_filename", "");
    min_lut   = param_tree.get("LUT.min_lut", true);
    
    
    input_string = param_tree.get_optional<std::string>("LUT.reuse_lut");
    if(input_string)    reuse_lut = bvec(*input_string);
    else reuse_lut = zeros_b(max_iter);
    
    
    
    
    // Set path objects and create folders if necessarry
    boost::system::error_code mkdir_error;
    
    trees_path = trees_dir;
    if(trees_path.is_relative()){
        trees_path = base_dir_path / trees_path;
    }
    create_directories( trees_path, mkdir_error );
    if(mkdir_error){
        it_error("LDPC_BER_Sim::LDPC_BER_Sim(): Could not create directory" << trees_path.string());
    }
    
}



void LDPC_BER_Sim_LUT::load(){
    
    if(codec_filename.length() == 0){ //design codec
        path parity_path = codes_path / (parity_filename + ".alist");
        if(!exists(parity_path)){
            it_error("Parity file" << parity_path.string() << " does not exist!");
        }
        
        // Load Parity check Matrix
        LDPC_Parity* H = new LDPC_Parity(parity_path.string());
        LDPC_Generator_Systematic* G = nullptr;
        
        // Load Generator if needed and present
        if(zero_codeword == false){
            path generator_path = codes_path / (parity_filename + ".gen.it");
            if(exists(generator_path)){
                G = new LDPC_Generator_Systematic();
                G->load(generator_path.string());
            }
            else{
                it_info_debug("No generator file found" << std::endl
                              << "Calculating Generator based on Parity Check Matrix ..." << std::endl);
                G = new LDPC_Generator_Systematic(H);
                
                // Save newly generated Generator and correspondingly permuted Parity
                if(save_permuted){
                    H->save_alist(parity_path.string());
                    it_file f(generator_path.string());
                    f<< Name("Fileversion") << 2;
                    f.close();
                    G->save(generator_path.string());
                    it_info_debug("done");
                }
            }
            
            encoder_set = true;
        }
        
        // Load/ Build tree structures
        LDPC_Code_LUT* C_ldpc = new LDPC_Code_LUT(H, G);
        codeword_length = C_ldpc->get_nvar();
        dataword_length = C_ldpc->get_ninfo();
        decoder_set = true;
        
        double sigma2_design = 0;
        if(design_thr)
            sigma2_design = sqr(*design_thr);
        else if(design_SNRdB)
            sigma2_design = pow10(-(*design_SNRdB)/10)/(2*C_ldpc->get_rate());
        else
            it_error("LDPC_BER_Sim_LUT::load(): No design noise threshold specified");
        

        if( tree_mode == "auto_bin_balanced" || tree_mode == "auto_bin_high" || tree_mode=="root_only" ){
            C_ldpc->design_luts("auto_bin_balanced", get_empirical_ensemble(*H) , min_lut,  sigma2_design , max_iter, reuse_lut, Nq_Cha, Nq_Msg);
        }
        else if (tree_mode == "file"){
            it_assert(trees_filename.length() > 0, "LDPC_BER_Sim_LUT::design_lut_codec(): Specify tree file name!");
            path tree_file_path = trees_path / trees_filename;
            it_assert( exists(tree_file_path), "LDPC_BER_Sim_LUT::design_lut_codec(): Tree file " << tree_file_path.string() << " could not be located");
            C_ldpc->design_luts("filename=" + tree_file_path.string(), get_empirical_ensemble(*H) ,  min_lut,  sigma2_design , max_iter, reuse_lut, Nq_Cha, Nq_Msg);
        }
        else{
            it_error("LDPC_BER_Sim_LUT::load(): tree_mode " << tree_mode << " unknown");
        }
        
        C_ldpc->set_exit_conditions(max_iter, parity_check_iter, parity_check_iter);
        C = C_ldpc;
        
        // Save codec
        if(rand_seed == save_codec){
            path results_path_subdir = results_path / gen_filename();
            boost::system::error_code mkdir_error;
            create_directories( results_path_subdir, mkdir_error );
            if(mkdir_error){
                it_error("LDPC_BER_Sim::load(): Could not create directory" << results_path_subdir.string());
            }
            C_ldpc->save_code((results_path_subdir / "lut_codec.it").string());
        }
    }// end if codec_filename>0
    else //load codec
    {
        path codec_path = codes_path / (codec_filename);
        if(!exists(codec_path)){
            it_error("Codec file" << codec_path.string() << " does not exist!");
        }
        
        LDPC_Generator_Systematic* G  = new LDPC_Generator_Systematic();
        LDPC_Code_LUT* C_ldpc = new LDPC_Code_LUT(codec_path.string(), G);
        C = C_ldpc;
        decoder_set = true;
        if (G->is_initialized() == true){
            encoder_set = true;
        }
        else{
            delete G;
            G = nullptr;
        }
        codeword_length = C_ldpc->get_nvar();
        dataword_length = C_ldpc->get_ninfo();
        max_iter = C_ldpc->get_nrof_iterations();
        Nq_Msg = C_ldpc->Nq_Msg;
        reuse_lut = C_ldpc->reuse_vec;
    }
    
    
    // Create results
    results = LDPC_BER_Sim_Results(codeword_length, codeword_length - dataword_length);
}


std::string LDPC_BER_Sim_LUT::gen_filename() const {
    it_assert(decoder_set, "LDPC_BER_Sim_LUT::gen_filename(): Decoder has not been set!");
    std::stringstream fn;
    fn  << results_prefix
    << "_N" << codeword_length
    << "_R" << C->get_rate()
    << "_maxIter" << max_iter
    << "_zcw" << (int) zero_codeword
    << "_frames" << Nframes;
    if(min_lut)
        fn << "_minLUT";
    else
        fn << "_LUT";
    fn <<  custom_name;
    return fn.str();
}
// ----------------------------------------------------------------------
// Other Functions
// ----------------------------------------------------------------------



