/*!
 * \file
 * \brief This is a simple statistical test wether the node updates are symmetric as they should be
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
    
    //==== Set Parameters
    
    int N_updates = 1e5;
    int Nq_fine = 5e3;
    int dv = 6;
    int dc = 32;
    int Nq_Msg = 3;
    int Nq_Cha = 4;
    double sig = .01;
    LUT_Tree var_tree = LUT_Tree("riiim/m//im/m///m//c/", LUT_Tree::VARTREE);
    LUT_Tree dec_tree = LUT_Tree("rim/m/m//im/m/m//c/", LUT_Tree::DECTREE);
    LUT_Tree chk_tree = LUT_Tree(dc-1, LUT_Tree::CHKTREE);
    std::cout << chk_tree.gen_template_string() << std::endl;
    
    
    //==== Get initial pmfs and design quantizers
    vec pmf_Cha, pmf_Msg;
    
    double delta = 2*(2/sqr(sig) + 4*2/sig)/Nq_fine;
    vec pmf_channel_fine = get_gaussian_pmf(2/sqr(sig), 2/sig, Nq_fine, delta);
    ivec Q_out;
    (void) quant_mi_sym(pmf_Cha,     Q_out, pmf_channel_fine, pow2i(Nq_Cha), true);
    (void) quant_mi_sym(pmf_Msg, Q_out, pmf_channel_fine, pow2i(Nq_Msg), true);
    
    var_tree.set_resolution(pow2i(Nq_Msg), pow2i(Nq_Msg), pow2i(Nq_Cha));
    var_tree.set_leaves(pmf_Msg, pmf_Cha);
    var_tree.update();
    
    dec_tree.set_resolution(pow2i(Nq_Msg), 2, pow2i(Nq_Cha));
    dec_tree.set_leaves(pmf_Msg, pmf_Cha);
    dec_tree.update();
    
    
    
    
    //==== Implementation
    
    for(int nn=0; nn<N_updates; nn++){
        int llr_in, llr_in_inv;
        int llr_out, llr_out_inv;
        ivec msg_in, msg_in_inv;
        ivec msg_out, msg_out_inv;

        msg_in = randi(dv ,0, pow2i(Nq_Msg)-1);
        msg_in_inv = pow2i(Nq_Msg)-1 - msg_in;
        
        llr_in = randi(0, pow2i(Nq_Cha)-1);
        llr_in_inv = pow2i(Nq_Cha)-1 - llr_in;
        
        std::deque<int> msg_q, msg_q_inv, msg_q_in, msg_q_in_inv;
        for(int ii=0; ii< dv; ii++){
            msg_q_in.push_back(msg_in(ii));
            msg_q_in_inv.push_back(msg_in_inv(ii));
        }
        msg_q = msg_q_in;
        msg_q_inv = msg_q_in_inv;
        msg_out     = var_tree.var_msg_update(msg_q, llr_in);
        msg_out_inv = var_tree.var_msg_update(msg_q_inv, llr_in_inv);
        
        msg_q = msg_q_in;
        msg_q_inv = msg_q_in_inv;
        llr_out     = dec_tree.dec_update(msg_q, llr_in);
        llr_out_inv = dec_tree.dec_update(msg_q_inv, llr_in_inv);
        
//        std::cout << "input msgs = " << msg_in << ", llr = " << llr_in << std::endl
//        << "out norm = " << msg_out << std::endl
//        << "out inve = " << msg_out_inv << "\n" << std::endl;
        
        if( sum(to_ivec((msg_out + msg_out_inv) ==  (pow2i(Nq_Msg)-1) )) !=  dv ){
            it_error("Non symmetric Variabe node update!");
	}
        
        if( llr_out + llr_out_inv != 1 ){
            it_error("Non symmetric Decision node update!");
	}
        
    }
    
    return EXIT_SUCCESS;
}
