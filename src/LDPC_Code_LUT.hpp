/*!
 * \file LDPC_Code_LUT.hpp
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


#ifndef LDPC_Code_LUT_hpp
#define LDPC_Code_LUT_hpp

#include "LDPC_DE.hpp"

#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iomanip>




namespace lut_ldpc{
    
    
/*!
 \brief LDPC decoder based on the LUT decoding
 
 LDPC decoder based on the LUT decoding and the min-LUT algorithm( cf. M. Meidlinger, A. Balatsoukas-Stimming, A. Burg, and G. Matz, “Quantized message passing for LDPC codes,” in Proc. 49th Asilomar Conf. Signals, Systems and Computers, Pacific Grove, CA, USA, Nov. 2015.
 )
 The following conventions are used:
    - At iteration ii, the quantizer for the LUT update has input resolution \c Nq_Msg(ii-1) and
      output resolution Nq(ii), where ii= 0,...,max_iters-1. For the first iteration(ii=0),
      the input resolution is implicityl given by the resolution of the channel LLR.
      Eventually, the last element Nq(max_iters-1) indicates the resolution of the output LLR.
      If we are not interested in iterative decoding-demodulation but only in the hard decoder
      output, this must be set to two.
    - If reuse_vec(ii+1) = true, then the LUT from iteration ii is reused at iteration ii+1 where ii=0,...,max_iters.
      For the first iteration (ii=0), a reuse is not possible and thus it must be that reuse_vec(0)=0. At the last iteration,
      it is pointless to perfor a variable node update because the newly calculated messages are not used for subsequent iteration.
      So at the last iteration (ii=max_iters-1), decoding is terminated after the check node update and the output LLR is calculated
      using decision node LUTs. For the design of the LUTs though, it is important that reuse_vec(max_iter-1)=0, so that the LUT for the last iteration
      is designed differently (namely as a decision LUT) rather that a message update LUT.
 */
class LDPC_Code_LUT : public Channel_Code {
    
friend class LDPC_BER_Sim_LUT;
    
private:
    /*!
     \brief
     How the messages for the initial decoding iteration are obtained.
     
     Options are "CONT" for messages from continuous inputs via quantization and
     "QCHA" for messages derived from the already quantized channel llrs"
     */
    enum initial_message_mode_e{
        CONT,
        QCHA,
        num_initial_message_modes
    }typedef initial_message_mode_t;
    //! How the messages for the initial decoding iteration are obtained
    initial_message_mode_t initial_message_mode;
    
public:
    //! Default constructor
    LDPC_Code_LUT();
    
    //! Constructor for Codec without LUTs
    LDPC_Code_LUT(const LDPC_Parity* const H,
                  LDPC_Generator* const G = 0,
                  bool perform_integrity_check = true);
    
    /*!
     \brief MinLUT Decoder Constructor, from a parity check matrix and optionally a
     generator.
     
     This constructor simply calls \c set_code().
     
     \param H The parity check matrix
     \param G A pointer to the optional generator object
     \param perform_integrity_check if true, then check that the parity and generator matrices are consistent
     */
    LDPC_Code_LUT(   const LDPC_Parity* const H,
                     const Array<Array<LUT_Tree>> var_trees_,
                     const bvec& reuse_vec_,
                     int Nq_Cha_,
                     const ivec& Nq_Msg_,
                     vec qb_Cha_,
                     vec qb_Msg_,
                     LDPC_Generator* const G = 0,
                     bool perform_integrity_check = true);
    /*!
     \brief LUT Decoder Constructor, from a parity check matrix and optionally a
     generator.
     
     This constructor simply calls \c set_code().
     
     \param H The parity check matrix
     \param G A pointer to the optional generator object
     \param perform_integrity_check if true, then check that the parity and generator matrices are consistent
     */
    LDPC_Code_LUT(   const LDPC_Parity* const H,
                  const Array<Array<LUT_Tree>> var_trees_,
                  const Array<Array<LUT_Tree>> chk_trees_,
                  const bvec& reuse_vec_,
                  int Nq_Cha_,
                  const ivec& Nq_Msg_,
                  vec qb_Cha_,
                  vec qb_Msg_,
                  LDPC_Generator* const G = 0,
                  bool perform_integrity_check = true);
    

    //! Constructor from filename to load pregenerated codec
    LDPC_Code_LUT(const std::string& filename, LDPC_Generator* const G_in=0):
        H_defined(false),
        G_defined(false),
        psc(true), pisc(false)
    {
        load_code(filename, G_in);
    };
    //! Destructor
    virtual ~LDPC_Code_LUT() {}
    
    /*!
     \brief Set the codec, from a parity check matrix and optionally a
     generator
     
     \param H The parity check matrix
     \param G A pointer to the optional generator object
     \param perform_integrity_check if true, then check that the parity and generator matrices are consistent
     */
    void set_code(const LDPC_Parity* const H, LDPC_Generator* const G = 0,
                  bool perform_integrity_check = true);

    
    //! Set the tree structure of the Min-LUT decoder
    void set_trees(const Array<Array<LUT_Tree>>& var_trees_,  bool perform_integrity_check=true);
    //! Set the tree structure of the Min-LUT decoder
    void set_trees(const Array<Array<LUT_Tree>>& var_trees_, const Array<Array<LUT_Tree>>& chk_trees_ ,  bool perform_integrity_check=true);
    
    //! Define the general LUT Tree structure (\c tree_method) and populate the structure with LUTs (\c lut_method)
    double design_luts(const std::string& tree_method,
                     const LDPC_Ensemble& ens,
                     bool min_lut,
                     double sigma2,
                     int max_iters,
                     const bvec& reuse_vec,
                     int Nq_Cha,
                     const ivec& Nq_Msg);
    
    
    /*!
     \brief Set the decoding loop exit conditions
     
     \param max_iters Maximum number of the decoding iterations
     \param syndr_check_each_iter If true, break the decoding loop as soon
     as valid codeword is found. Recommended value: \c true.
     \param syndr_check_at_start If true, perform a syndrome check before
     entering the decoding loop. If LLRin corresponds to a valid codeword,
     set LLRout = LLRin. Recommended value: \c false.
     
     \note The default values set in the class constructor are: "50",
     "true" and "false", respectively.
     */
    void set_exit_conditions(int max_iters,
                             bool syndr_check_each_iter = true,
                             bool syndr_check_at_start = false);
    
    
    
    // ------------ Encoding  ---------------------
    
    /*!
     \brief Encode codeword
     
     This is a wrapper functions which calls a proper implementation from
     the \c LDPC_Generator object.
     
     \param input Vector of \c ncheck input bits
     \param output Vector of \c nvar output bits
     */
    virtual void encode(const bvec &input, bvec &output); // Implemented in base class
    //! \brief Encode codeword
    virtual bvec encode(const bvec &input); // Implemented in base class
    
    
    // ------------ Decoding  ---------------------
    
    //! Inherited from the base class - not implemented here
    virtual void decode(const bvec &, bvec &) {
        it_error("LDPC_Code_MinLUT::decode(): Hard input decoding not implemented");
    }
    //! Inherited from the base class - not implemented here
    virtual bvec decode(const bvec &) {
        it_error("LDPC_Code_MinLUT::decode(): Hard input decoding not implemented");
        return bvec();
    }
    
    //! This function outputs systematic bits of the decoded codeword
    virtual void decode(const vec &llr_in, bvec &syst_bits);
    //! This function outputs systematic bits of the decoded codeword
    virtual bvec decode(const vec &llr_in);
    
    //! Inherited from the base class - not implemented here
    void decode_soft_out(const vec &llr_in, vec &llr_out);
    //! Inherited from the base class - not implemented here
    vec decode_soft_out(const vec &llr_in);
    
    //! Implements the lut decoding algorithm
    int lut_decode(const ivec LLRin_cha, const ivec LLRin_msg,  bvec& LLRout);
    
    //! Performs a syndrom check based the state of the check to variable node messages. For this, it needs to know the resolution of the messages to decide which ones are negative. The binary vector \c b returns the unanimous message decision.
    bool syndrome_check(int Nq_Msg, bvec& b ) const;
    //! Syndrom Check based on binary vector \c b
    bool syndrome_check(const bvec &b) const;
    

    
    
    // ------------ Basic information gathering functions ------
    
    //! Get the coderate
    double get_rate() const {
        return (1.0 - static_cast<double>(nchk_lin_indep) / nvar);
    }
    
    //! Get the number of variable nodes
    int get_nvar() const { return nvar; }
    
    //! Get the number of check nodes
    int get_nchk() const { return nchk; }
    //! Get the number of check nodes
    int get_ncheck() const { return nchk; }
    
    //! Get the number of linearly independent check nodes
    int get_nchk_lin_indep() const { return nchk_lin_indep; }
    //! Get the number of linearly independent  check nodes
    int get_ncheck_lin_indep() const { return nchk_lin_indep; }
    
    //! Get the number of information bits per codeword
    int get_ninfo() const { return nvar - nchk_lin_indep; }
    
    //! Get the maximum number of LDPC decoding iterations
    int get_nrof_iterations() const { return max_iters; }
    
    //! Set decoder output verbosity
    void set_output_verbosity(int verb_level){ this->output_verbosity = verb_level;};
    //! Set initial message mode
    void set_initial_message_mode(initial_message_mode_t m){ this->initial_message_mode = m;};
    
    
    inline void chk_update_minsum(int node_idx, int edge_idx, int iter);
    inline void chk_update_lut(int node_idx, int edge_idx, int iter);
    inline void var_update_lut(int node_idx, int edge_idx, int iter, int llr);
    inline int dec_update_lut(int node_idx, int edge_idx, int iter, int llr);
    
    //! Save LUT Codec to file
    void save_code(const std::string& filename) const;
    //! Load LUT Codec from file
    void load_code(const std::string& filename, LDPC_Generator* const G = 0);
    
    friend  std::ostream &operator<<(std::ostream &os, const LDPC_Code_LUT &C);
    
protected:

    void decoder_parameterization(const LDPC_Parity* const H);
    void integrity_check();
    void setup_decoder();
    
private:
    bool H_defined;  //!< true if parity check matrix is defined
    bool G_defined;  //!< true if generator is defined
    bool LUTs_defined;  //! true if the LUT_Tree arrays are defined
    bool minLUT; //!< true if the checknode update is adhering to the min-sum algorithm
    int nvar;   //!< Number of variable nodes
    int nchk;   //!< Number of check nodes
    int nchk_lin_indep;   //!< Number of linearly independent check nodes
    LDPC_Generator *G;  //!< Generator object pointer
    
    // decoder parameters
    int max_iters;  //!< Maximum number of iterations
    bool psc;   //!< check syndrom after each iteration
    bool pisc;   //!< check syndrom before first iteration

    
    ivec dv_vec, dc_vec;
    int num_edges;
    ivec msgs;
    ivec cn_msg_idx;
    Array<ivec> chk_equ_idx;
    
    //! LUT Tree array for variable node updates
    Array<Array<LUT_Tree>> var_trees;
    //! LUT Tree array for check node updates
    Array<Array<LUT_Tree>> chk_trees;
    
    //! Vector of length \c max_iters, indicating which element from \c var_trees should be used
    ivec var_tree_idx_iter;
    //! Vector of length \c max_iters, indicating which element from \c chk_trees should be used
    ivec chk_tree_idx_iter;
    //! Size of the first dimension of \c var_trees
    int num_var_trees_iter;
    //! Size of the first dimension of \c chk_trees
    int num_chk_trees_iter;
    //! Reuse vector
    bvec reuse_vec;

    
    //! Resolution of channel input
    int Nq_Cha;
    //! Resolution of messages for each iteration
    ivec Nq_Msg;
    
    //! Quantization boundaries to transform continous input LLRs to a discrete internal LLR representation
    vec qb_Cha, qb_Msg;
    
    /*!
     \brief
     Information optimal quantizer for mapping discrete messages from Nq_Cha resolution to Nq_Msg resolution
     
     This vector is not used internally for decoding but is rather used when exporting the decoder to VHDL
     */
    ivec Nq_Cha_2_Nq_Msg_map;
    
    //! If != 0, pairs of decoder inputs and outputs are printed to the standard output.
    int output_verbosity;
    
    
    
    //! Vector of length nvar indicating which tree to use for the var updates
    ivec var_tree_idx_degree;
    //! Vector of length nchk indicating which tree to use for the chk updates
    ivec chk_tree_idx_degree;
    //! Size of the second dimension of \c var_trees
    int num_var_trees_degree;
    //! Size of the second dimension of \c chk_trees
    int num_chk_trees_degree;
    
    //! Maximum check node degree that the class can handle
    static const int max_cnd = 200;

 
};
    
//! 
std::ostream &operator<<(std::ostream &os, const LDPC_Code_LUT &C);
}

#endif
