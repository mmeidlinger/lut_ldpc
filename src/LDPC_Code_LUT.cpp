/*!
 * \file LDPC_Code_LUT.cpp
 * \brief Implementation of LUT based LDPC encoding and decoding for
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


#include "LDPC_Code_LUT.hpp"

using namespace itpp;
using namespace lut_ldpc;

static const int LUT_LDPC_binary_file_version = 1;

// ----------------------------------------------------------------------
// LDPC_Code_LUT
// ----------------------------------------------------------------------
LDPC_Code_LUT::LDPC_Code_LUT():
    initial_message_mode(CONT),
    H_defined(false),
    G_defined(false),
    LUTs_defined(false),
    max_iters(0),
    psc(true),
    pisc(false),
    output_verbosity(0)
{}

LDPC_Code_LUT::LDPC_Code_LUT(const LDPC_Parity* const H,
                             LDPC_Generator* const G,
                             bool perform_integrity_check):
initial_message_mode(CONT), H_defined(false), G_defined(false), LUTs_defined(false), output_verbosity(0)
{
    
    set_code(H, G, perform_integrity_check);
}
    
LDPC_Code_LUT::LDPC_Code_LUT(const LDPC_Parity* const H_,
                     const Array<Array<LUT_Tree>> var_trees_,
                     const bvec& reuse_vec_,
                     int Nq_Cha_,
                     const ivec& Nq_Msg_,
                     vec qb_Cha_,
                     vec qb_Msg_,
                     LDPC_Generator* const G_,
                     bool perform_integrity_check):
initial_message_mode(CONT), H_defined(false), G_defined(false), LUTs_defined(false), output_verbosity(0)
{
    
    
    
    set_code(H_, G_, perform_integrity_check);
    
    reuse_vec = reuse_vec_;
    max_iters = length(reuse_vec);
    Nq_Cha = Nq_Cha_;
    Nq_Msg = Nq_Msg_;
    qb_Cha = qb_Cha_;
    qb_Msg = qb_Msg_;
    
    set_trees(var_trees_, perform_integrity_check);
    minLUT = true;
    LUTs_defined = true;

}
    
    
LDPC_Code_LUT::LDPC_Code_LUT(const LDPC_Parity* const H_,
                             const Array<Array<LUT_Tree>> var_trees_,
                             const Array<Array<LUT_Tree>> chk_trees_,
                             const bvec& reuse_vec_,
                             int Nq_Cha_,
                             const ivec& Nq_Msg_,
                             vec qb_Cha_,
                             vec qb_Msg_,
                             LDPC_Generator* const G_,
                             bool perform_integrity_check){
    
    
    
    set_code(H_, G_, perform_integrity_check);
    
    reuse_vec = reuse_vec_;
    max_iters = length(reuse_vec);
    Nq_Cha = Nq_Cha_;
    Nq_Msg = Nq_Msg_;
    qb_Cha = qb_Cha_;
    qb_Msg = qb_Msg_;
    
    set_trees(var_trees_, chk_trees_,  perform_integrity_check);
    minLUT = false;
    LUTs_defined = true;
    output_verbosity = 0;
    initial_message_mode = CONT;
    
}

void LDPC_Code_LUT::set_trees(const Array<Array<LUT_Tree>>& var_trees_, const Array<Array<LUT_Tree>>& chk_trees_ ,  bool perform_integrity_check){
    
    it_assert(reuse_vec(0) == bin(0) && reuse_vec(length(reuse_vec)-1) == bin(0), "LDPC_Code_LUT::set_trees(): First and last iteration are excempt from tree reuse");
    
    // Set var trees
    var_tree_idx_iter = cumsum(to_ivec(reuse_vec==bin(0)))-1;
    num_var_trees_iter = sum(to_ivec(reuse_vec == bin(0)));
    
    var_tree_idx_degree = zeros_i(nvar);
    num_var_trees_degree = var_trees_(0).size();
    
    ivec var_tree_degrees = zeros_i(num_var_trees_degree);
    
    for(int dd=0; dd<num_var_trees_degree; dd++){
        var_tree_degrees(dd) = var_trees_(0)(dd).get_num_leaves();
    }
    for(int nn=0; nn< nvar; nn++){
        ivec idx = find(var_tree_degrees == dv_vec(nn));
        var_tree_idx_degree(nn) = idx(0);
    }
    this->var_trees = var_trees_;
    
    // Set chk trees
    if(chk_trees_.size()>0){
        chk_tree_idx_iter = cumsum(to_ivec(reuse_vec==bin(0)))-1;
        num_chk_trees_iter = sum(to_ivec(reuse_vec == bin(0)));
        
        chk_tree_idx_degree = zeros_i(nchk);
        num_chk_trees_degree = chk_trees_(0).size();
        
        ivec chk_tree_degrees = zeros_i(num_chk_trees_degree);
        
        for(int dd=0; dd<num_chk_trees_degree; dd++){
            chk_tree_degrees(dd) = chk_trees_(0)(dd).get_num_leaves()+1;
        }
        for(int nn=0; nn< nchk; nn++){
            ivec idx = find(chk_tree_degrees == dc_vec(nn));
            chk_tree_idx_degree(nn) = idx(0);
        }
        this->chk_trees = chk_trees_;
    }
    
    // Check wether the trees
    //  * match the Parity check matrix (active degrees)
    //  * the tree resolutions match the vector of resolutions and are consistent (input resolution of ii+1 = output resolution of ii)
    //  *
    if(perform_integrity_check){
        
    }
}
 
void LDPC_Code_LUT::set_trees(const Array<Array<LUT_Tree>>& var_trees_,  bool perform_integrity_check){
    Array<Array<LUT_Tree>> chk_trees_dummy(0);
    set_trees(var_trees_, chk_trees_dummy,  perform_integrity_check);
}
    
void LDPC_Code_LUT::set_exit_conditions(int max_iters_in,
                                    bool syndr_check_each_iter,
                                    bool syndr_check_at_start)
{
    it_assert(max_iters >= 0, "LDPC_Code_LUT::set_nrof_iterations(): Maximum "
              "number of iterations can not be negative");
    max_iters = max_iters_in;
    psc = syndr_check_each_iter;
    pisc = syndr_check_at_start;
}

    
void LDPC_Code_LUT::encode(const bvec &input, bvec &output)
{
    it_assert(G_defined, "LDPC_Code_LUT::encode(): LDPC Generator is required "
              "for encoding");
    G->encode(input, output);
    it_assert_debug(syndrome_check(output), "LDPC_Code::encode(): Syndrome "
                    "check failed");
}

bvec LDPC_Code_LUT::encode(const bvec &input)
{
    bvec result;
    encode(input, result);
    return result;
}

void LDPC_Code_LUT::decode(const vec &llr_in, bvec &syst_bits)
{
    // Get quantized channel llrs
    ivec llr_in_cha = quant_nonlin(llr_in, qb_Cha);
    
    // Get quantized messages for the initial decoding iteration
    ivec llr_in_msg;
    switch (this->initial_message_mode) {
        case CONT:
            llr_in_msg = quant_nonlin(llr_in, qb_Msg);
            break;
        case QCHA:
            llr_in_msg = Nq_Cha_2_Nq_Msg_map(llr_in_cha);
            break;
        default:
            it_error("LDPC_Code_LUT::decode(): Initial message mode undefined!");
            break;
    }
    
    bvec llr_out;
    lut_decode(llr_in_cha, llr_in_msg,  llr_out);
    syst_bits = llr_out.left(nvar - nchk_lin_indep);
    
    // Print input and output pairs
    if(output_verbosity>0){
        std::cout << "Stimuli Pair (Quantized channel LLR decoder inputs in hex format and decoder output in binary format): " << std::endl;
        for(int ii=0; ii<llr_in_cha.length(); ii++){
            std::cout << std::setfill('0') << std::setw(8) << std::uppercase << std::hex << llr_in_cha(ii) << "  ";
        }
        std::cout << std::endl;
        for(int ii=0; ii<llr_in_cha.length(); ii++){
            std::cout << llr_out(ii) << "  ";
        }
        std::cout << std::endl << std::endl;
    }
}

bvec LDPC_Code_LUT::decode(const vec &llr_in)
{
    bvec syst_bits;
    decode(llr_in, syst_bits);
    return syst_bits;
}

void LDPC_Code_LUT::decode_soft_out(const vec &llr_in, vec &llr_out)
{
    it_error("LDPC_Code_LUT::decode_soft_out(): Soft output decoding not implemented possible for LUT-type decode");
}

vec LDPC_Code_LUT::decode_soft_out(const vec &llr_in)
{
    it_error("LDPC_Code_LUT::decode_soft_out(): Soft output decoding not implemented possible for LUT-type decode");
    return vec(0);
}
    
int LDPC_Code_LUT::lut_decode(const ivec LLRin_cha, const ivec LLRin_msg,  bvec& LLRout){
    // Note the IT++ convention that a sure zero corresponds to
    // LLR=+infinity
    it_assert(H_defined, "LDPC_Code_LUT::lut_decode(): Parity check matrix not "
              "defined");
    it_assert(LUTs_defined, "LDPC_Code_LUT::lut_decode(): LUTs not defined");
    it_assert((LLRin_cha.size() == nvar) &&
              (LLRin_msg.size() == nvar) &&
              (dv_vec.size() == nvar) &&
              (dc_vec.size() == nchk), "LDPC_Code_LUT::lut_decode(): Wrong "
              "input dimensions");
    
    // Todo: Check ranges of LLRin_cha and LLRin_msg against resolution
    
    
    // Initially set the output LLR according to the sign of the input llr
    LLRout = (LLRin_cha < Nq_Cha/2);
    
    if (pisc && syndrome_check(LLRout)) {
        return 0;
    }
    
    int edge_idx = 0;
    
    // initial step (setting msgs to zero is not necessary due to this step)
    for (int vv = 0; vv < nvar; vv++) {
        for (int ee = 0; ee < dv_vec(vv); ee++) {
            msgs(edge_idx) = LLRin_msg(vv);
            edge_idx++;
        }
    }
    
    // Print initial VN-to-CN messages
    if(output_verbosity>1){
        std::cout << "Initial VN-to-CN messages: " << std::endl;
        for (int ee = 0; ee < num_edges; ee++) {
            std::cout << std::setfill('0') << std::setw(8) << std::uppercase << std::hex << msgs(ee) << "  ";
        }
        std::cout << std::endl;
    }

    int ii;
    for(ii=0; ii<max_iters; ii++){
        // Step 1: CN update
        edge_idx = 0;
        for(int cc=0; cc<nchk; cc++){
            if(minLUT)
                chk_update_minsum(cc,edge_idx, ii);
            else
                chk_update_lut(cc, edge_idx, ii);
            edge_idx += dc_vec(cc);
        }
        if(output_verbosity>2){ // Print CN-to-VN messages
            std::cout << "CN-to-VN messages after CN update at iteration " << ii << ":" << std::endl;
            for (int ee = 0; ee < num_edges; ee++) {
                std::cout << std::setfill('0') << std::setw(8) << std::uppercase << std::hex << msgs(ee) << "  ";
            }
            std::cout << std::endl;
        }
        
        // Step 2: VN update
        if(ii!= max_iters-1){
            edge_idx = 0;
            for(int vv=0; vv<nvar; vv++){
                var_update_lut(vv,edge_idx, ii, LLRin_cha(vv));
                edge_idx += dv_vec(vv);
            }
            //parity check based on unanimous message decision. This
            if (psc && syndrome_check(Nq_Msg(ii+1),LLRout)) {
                return ii+1;
            }
        }
        if(output_verbosity>1){ // Print VN-to-CN messages
            std::cout << "VN-to-CN messages after VN update at iteration " << ii << ":" << std::endl;
            for (int ee = 0; ee < num_edges; ee++) {
                std::cout << std::setfill('0') << std::setw(8) << std::uppercase << std::hex << msgs(ee) << "  ";
            }
            std::cout << std::endl;
        }
    }
    // Decision node step
    edge_idx = 0;
    for(int vv=0; vv<nvar; vv++){
        LLRout(vv) = (dec_update_lut(vv,edge_idx, max_iters-1, LLRin_cha(vv)) < 1 );
        edge_idx += dv_vec(vv);
    }
    
    if (syndrome_check(LLRout))
        return max_iters;
    else
        return -max_iters;
    
    

}
 
inline void LDPC_Code_LUT::chk_update_minsum(int node_idx, int edge_idx, int iter){
    int min1,min2, tmp;
    int min_idx = 0;
    int sign_prod = 0; /*0 indicates a positive sign*/
    int sign_msg;
    
    int nz = Nq_Msg(iter)/2;
    
    min1 = nz;
    min2 = nz;
    
    /*Find Minimum and second minimum and calculate the product of all signs*/
    for(int cc=0; cc<dc_vec(node_idx); cc++){
        if(msgs[cn_msg_idx[edge_idx+cc]]< nz){
            sign_prod ^= 1;
            tmp =  nz-1 - msgs[cn_msg_idx[edge_idx+cc]];
        }
        else{
            tmp =  msgs[cn_msg_idx[edge_idx+cc]] - nz;
        }
        if(tmp < min1){
            min2    = min1;
            min1    = tmp;
            min_idx = cc;
        }
        else if(tmp < min2)
            min2 = tmp;
    }
    /*Assign output messages*/
    for(int cc=0; cc<dc_vec(node_idx); cc++){
        /* The magnitude of the output is either min1 or min2 */
        if(cc == min_idx)
            tmp = min2;
        else
            tmp = min1;
        /* Invert the impact of the own sign */
        if(msgs[cn_msg_idx[edge_idx+cc]]< nz )
            sign_msg = sign_prod ^ 1;
        else
            sign_msg = sign_prod;
        /* Assign output */
        if(sign_msg)
            msgs[cn_msg_idx[edge_idx+cc]] = nz-1 - tmp;
        else
            msgs[cn_msg_idx[edge_idx+cc]] = nz + tmp;
        
    }
}
    
inline void LDPC_Code_LUT::var_update_lut(int node_idx, int edge_idx, int iter, int llr){
    std::deque<int> msgs_in;
    for(int ii=0; ii< dv_vec(node_idx); ii++){
        msgs_in.push_back(msgs(edge_idx+ii));
    }
    ivec msgs_out = var_trees(var_tree_idx_iter(iter))(var_tree_idx_degree(node_idx)).var_msg_update(msgs_in, llr);
    for(int ii=0; ii< dv_vec(node_idx); ii++){
        msgs(edge_idx+ii) = msgs_out(ii);
    }
    return;
}

inline void LDPC_Code_LUT::chk_update_lut(int node_idx, int edge_idx, int iter){
    std::deque<int> msgs_in;
    for(int ii=0; ii< dc_vec(node_idx); ii++){
        msgs_in.push_back(msgs(cn_msg_idx[edge_idx+ii]));
    }
    ivec msgs_out = chk_trees(chk_tree_idx_iter(iter))(chk_tree_idx_degree(node_idx)).chk_msg_update(msgs_in);
    for(int ii=0; ii< dc_vec(node_idx); ii++){
        msgs(cn_msg_idx[edge_idx+ii]) = msgs_out(ii);
    }
    return;
}

inline int LDPC_Code_LUT::dec_update_lut(int node_idx, int edge_idx, int iter, int llr){
    std::deque<int> msgs_in;
    for(int ii=0; ii< dv_vec(node_idx); ii++){
        msgs_in.push_back(msgs(edge_idx+ii));
    }
    return  var_trees(var_tree_idx_iter(iter))(var_tree_idx_degree(node_idx)).dec_update(msgs_in, llr);
}


bool LDPC_Code_LUT::syndrome_check(int Nq_Msg, bvec& b) const {
    int edge_idx = 0;
    int nz = Nq_Msg/2;
    it_assert(b.size() == nvar, "LDPC_Code_LUT::syndrome_check(): Vector must be of size nvar!");
    
    for(int vv=0; vv<nvar; vv++){
        bool bit = ( msgs(edge_idx)< nz);
        for(int ee=1; ee<dv_vec(vv); ee++){
            bool tmp = (msgs(edge_idx+ee)< nz);
            if (bit != tmp )    return false;
        }
        edge_idx += dv_vec(vv);
        b(vv) = bit;
    }
    return syndrome_check(b);
}
    

bool LDPC_Code_LUT::syndrome_check(const bvec &b) const {
    // a sure zero corresponds to a  LLR +infinity
    for(int cc = 0; cc<nchk; cc++){
        int synd = 0;
        for (int ii = 0; ii < dc_vec(cc); ii++) {
            if (b(chk_equ_idx(cc)(ii))) {
                synd++;
            }
        }
        if ((synd&1) == 1) {
            return false;  // codeword is invalid
        }
    }
    return true;
}

void LDPC_Code_LUT::set_code(const LDPC_Parity* const H,
                         LDPC_Generator* const G_in,
                         bool perform_integrity_check)
{
    decoder_parameterization(H);
    G = G_in;
    if (G != 0) {
        G_defined = true;
        if (perform_integrity_check) {
            integrity_check();
        } else {
            it_info_debug("LDPC_Code::set_code(): integrity check was not performed");
        }
    }
}
// Private methods
    
void LDPC_Code_LUT::decoder_parameterization(const LDPC_Parity* const Hmat)
{
    // copy basic parameters
    nvar = Hmat->nvar;
    nchk = Hmat->ncheck;
    if(nvar < 1e5){
        nchk_lin_indep = GF2mat(full(Hmat->get_H())).row_rank();
    }
    else{
        it_info("LDPC_Code_LUT::decoder_parameterization(): Code too large to determine rank. Setting nchk_lin_indep = nchk");
        nchk_lin_indep = nchk;
    }
    
    dv_vec = Hmat->sumX1;
    dc_vec  = Hmat->sumX2;
    num_edges = sum(dv_vec);

    
    // Allocate idx vectors/ array
    cn_msg_idx= ivec(num_edges);
    msgs = ivec(num_edges);
    chk_equ_idx = Array<ivec>(nchk);
    
    
    // Generate check node indices
    Array<ivec> tmp = Array<ivec>(nchk);
    int c=0;
    for(int vv=0;vv<nvar;vv++){
        ivec idx = Hmat->get_col(vv).get_nz_indices();
        sort(idx);
        for(int ii=0;ii< idx.length(); ii++){
            tmp(idx(ii)) = concat(tmp(idx(ii)), c);
            c++;
        }
    }
    c=0;
    for(int cc=0;cc<nchk;cc++){
        cn_msg_idx.set_subvector(c, tmp(cc));
        c += tmp(cc).length();
    }
    it_assert(c == num_edges, "LDPC_Code_LUT::decoder_parameterization(): Dimension mismatch");
    
    // Save the index for the parity check equations
    for(int cc=0; cc<nchk; cc++){
        ivec idx = Hmat->get_row(cc).get_nz_indices();
        sort(idx);
        chk_equ_idx(cc) = idx;
    }
    

    
    H_defined = true;
    return;
}





void LDPC_Code_LUT::integrity_check()
{
    if (G_defined) {
        it_info_debug("LDPC_Code_LUT::integrity_check(): Checking integrity of "
                      "the LDPC_Parity and LDPC_Generator data");
        bvec bv(nvar - nchk_lin_indep), cw;
        bv.clear();
        bv(0) = 1;
        for (int i = 0; i < nvar - nchk_lin_indep; i++) {
            G->encode(bv, cw);
            it_assert(syndrome_check(cw),
                      "LDPC_Code_LUT::integrity_check(): Syndrome check failed");
            bv.shift_right(bv(nvar - nchk_lin_indep - 1));
        }
    }
    else {
        it_info_debug("LDPC_Code_LUT::integrity_check(): No generator defined "
                      "- no check performed");
    }
}

void LDPC_Code_LUT::load_code(const std::string& filename, LDPC_Generator* const G_in)
{
    it_info_debug("LDPC_Code_LUT::load_code(): Loading LDPC LUT codec from "
                  << filename);
    
    it_ifile f(filename);
    int ver;
    std::string lut_string;
    std::stringstream lut_stream;
    
    f >> Name("Fileversion") >> ver;
    it_assert(ver == LUT_LDPC_binary_file_version, "LDPC_Code_LUT::load_code(): Unsupported file format");
    f >> Name("H_defined") >> H_defined;
    f >> Name("G_defined") >> G_defined;
    f >> Name("LUTs_defined") >> LUTs_defined;
    f >> Name("nvar") >> nvar;
    f >> Name("nchk") >> nchk;
    f >> Name("nchk_lin_indep") >> nchk_lin_indep;
    f >> Name("dv_vec") >> dv_vec;
    f >> Name("dc_vec") >> dc_vec;
    f >> Name("chk_equ_idx") >> chk_equ_idx;
    f >> Name("cn_msg_idx") >> cn_msg_idx;
    f >> Name("max_iters") >> max_iters;
    
    f >> Name("Nq_Cha") >> Nq_Cha;
    f >> Name("Nq_Msg") >> Nq_Msg;
    f >> Name("Nq_Cha_2_Nq_Msg_map") >> Nq_Cha_2_Nq_Msg_map;
    f >> Name("qb_Cha") >> qb_Cha;
    f >> Name("qb_Msg") >> qb_Msg;
    f >> Name("reuse_vec") >> reuse_vec;
    f >> Name("minLUT") >> minLUT;
    f >> Name("output_verbosity") >> output_verbosity;
    
    f >> Name("var_tree_string") >> lut_string;
    lut_stream << lut_string;
    lut_stream >> var_trees;
    
    lut_stream.clear();
    lut_string.clear();
    lut_stream.str(std::string());
    
    
    f >> Name("chk_tree_string") >> lut_string;
    lut_stream << lut_string;
    lut_stream >> chk_trees;
    
    f.close();
    
    // Set trees
    set_trees(var_trees, chk_trees, true);
    

    // Alocate message memory
    num_edges = sum(dv_vec);
    msgs = ivec(num_edges);
    
    // load generator data
    if (G_defined) {
        it_assert(G_in != 0, "LDPC_Code_LUT::load_code(): Generator object is "
                  "missing. Can not load the generator data from a file.");
        G = G_in;
        G->load(filename);
    }
    else {
        G = 0;
        it_info_debug("LDPC_Code_LUT::load_code(): Generator data not loaded. "
                      "Generator object will not be used.");
    }
    
    it_info_debug("LDPC_Code_LUT::load_code(): Successfully loaded LDPC LUT codec "
                  "from " << filename);
    

}

void LDPC_Code_LUT::save_code(const std::string& filename) const
{
    it_assert(H_defined, "LDPC_Code_LUT::save_to_file(): There is no parity "
              "check matrix");
    it_info_debug("LDPC_Code_LUT::save_to_file(): Saving LDPC codec to "
                  << filename);
    
    std::string lut_string;
    std::stringstream lut_stream;
    
    it_file f;
    f.open(filename, true);
    f << Name("Fileversion") << LUT_LDPC_binary_file_version;
    f << Name("H_defined") << H_defined;
    f << Name("G_defined") << G_defined;
    f << Name("LUTs_defined") << LUTs_defined;
    f << Name("nvar") << nvar;
    f << Name("nchk") << nchk;
    f << Name("nchk_lin_indep") << nchk_lin_indep;
    f << Name("dv_vec") << dv_vec;
    f << Name("dc_vec") << dc_vec;
    f << Name("chk_equ_idx") << chk_equ_idx;
    f << Name("cn_msg_idx") << cn_msg_idx;
    
    f << Name("Nq_Cha") << Nq_Cha;
    f << Name("Nq_Msg") << Nq_Msg;
    f << Name("Nq_Cha_2_Nq_Msg_map") << Nq_Cha_2_Nq_Msg_map;
    f << Name("qb_Cha") << qb_Cha;
    f << Name("qb_Msg") << qb_Msg;
    f << Name("reuse_vec") << reuse_vec;
    f << Name("minLUT") << minLUT;
    f << Name("output_verbosity") << output_verbosity;
    f << Name("max_iters") << max_iters;
    
    lut_stream << var_trees;
    f << Name("var_tree_string") << lut_stream.str();
    
    lut_stream.clear();
    lut_stream.str(std::string());
    
    lut_stream << chk_trees;
    f << Name("chk_tree_string") << lut_stream.str();
    
    f.close();
    
    // save generator data;
    if (G_defined)
        G->save(filename);
    else
        it_info_debug("LDPC_Code_LUT::save_code(): Missing generator object - "
                      "generator data not saved");
    
    it_info_debug("LDPC_Code_LUT::save_code(): Successfully saved LDPC codec to "
                  << filename);
}

double LDPC_Code_LUT::design_luts(const std::string& tree_method,
                                const LDPC_Ensemble& ens,
                                bool min_lut,
                                double sigma2,
                                int max_iters,
                                const bvec& reuse_vec,
                                int Nq_Cha,
                                const ivec& Nq_Msg){

    this->minLUT = min_lut;
    this->max_iters = max_iters;
    this->reuse_vec = reuse_vec;
    this->Nq_Cha = Nq_Cha;
    this->Nq_Msg = Nq_Msg;

    double sig;
    
    Array<Array<LUT_Tree> > var_luts;
    Array<Array<LUT_Tree> > chk_luts;
    get_lut_tree_templates(tree_method, ens, Nq_Msg, Nq_Cha, min_lut, var_luts, chk_luts );
    
    
    // Generate LUTs and set object properties

    // Design LUTs with DE
    LDPC_DE_LUT de(ens, Nq_Cha, Nq_Msg, max_iters, var_luts, chk_luts, reuse_vec);

    sig = std::sqrt(sigma2);

    
    de.get_quant_bound(sig, qb_Cha, qb_Msg);
    de.get_lut_trees(var_luts, chk_luts, sig);

    // Set codec properties
    set_trees(var_luts, chk_luts);
    
    // Calculate Nq_Cha_2_Nq_Msg_map
    const double LLR_max_mag = 25.0;
    double delta = 2*LLR_max_mag/Nq_Cha;
    vec p_cha = get_gaussian_pmf(2/sigma2, 2/sig, Nq_Cha, 0);
    vec pmf_channel= get_gaussian_pmf(2/sqr(sig), 2/sig, Nq_Cha, delta);
    vec p_msg;
    (void) quant_mi_sym(p_msg, Nq_Cha_2_Nq_Msg_map, pmf_channel, Nq_Msg(0), true);
  
    LUTs_defined = true;
    return sig;
    
}
    
// ----------------------------------------------------------------------
// Related functions
// ----------------------------------------------------------------------

std::ostream& lut_ldpc::operator<<(std::ostream &os, const LDPC_Code_LUT &C)
{
    ivec rdeg = zeros_i(max(C.dc_vec) + 1);
    for (int i = 0; i < C.nchk; i++)     {
        rdeg(C.dc_vec(i))++;
    }
    
    ivec cdeg = zeros_i(max(C.dv_vec) + 1);
    for (int j = 0; j < C.nvar; j++)     {
        cdeg(C.dv_vec(j))++;
    }
    
    os << "--- LDPC codec ----------------------------------\n"
    << "Nvar : " << C.get_nvar() << "\n"
    << "Ncheck : " << C.get_nchk() << "\n"
    << "Rate : " << C.get_rate() << "\n"
    << "Column degrees (node perspective): " << cdeg << "\n"
    << "Row degrees (node perspective): " << rdeg << "\n"
    << "-------------------------------------------------\n"
    << "Decoder parameters:\n"
    << " - max. iterations : " << C.max_iters << "\n"
    << " - syndrome check at each iteration : " << C.psc << "\n"
    << " - syndrome check at start : " << C.pisc << "\n"
    << "-------------------------------------------------\n"
    << "Quantizer Settings: (TODO)\n";
    return os;
}




