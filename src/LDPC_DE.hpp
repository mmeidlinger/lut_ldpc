/*!
 * \file
 * \brief
 * \author Michael Meidlinger
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2016 Michael Meidlinger - All Rights Reserved
 *
 */



#ifndef LDPC_DE_hpp
#define LDPC_DE_hpp

#include "LUT_Tree.hpp"
#include "LDPC_Ensemble.hpp"



#include <iostream>
#include <itpp/itbase.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iomanip>

//! Git version of program
extern const char *gitversion;
using namespace itpp;

namespace lut_ldpc{


    

class LDPC_DE
{
public:
    enum {
        ARI,
        GEO
    };
    
    //! Returns the number of iterations if \c thr is feasible for \c ens, otherwise returns -1
    virtual int evolve(double thr) = 0;
    virtual int evolve(double thr, bool var, bool chk, mat& P, vec& p) = 0;
    //! This function is used to set and update the ensemble of the DE object
 //   virtual void setup(const LDPC_Ensemble &ens) = 0;
    
    
    //! Performs a bisection search and returns the threshold of \c ens
    virtual int bisec_search(double & sig);
    void set_bisec_window(double tmin, double tmax);
    
    //! Set the ensemble
    virtual void set_ensemble(const LDPC_Ensemble& ens) = 0;
    
    //! Set Exit conditions
    virtual void set_exit_conditions(int maxiter_de, int maxiter_bisec, int max_ni_de_iters, double Pe_max, double thr_prec){
        this->maxiter_de = maxiter_de;
        this->maxiter_bisec = maxiter_bisec;
        this->max_ni_de_iters = max_ni_de_iters;
        this->Pe_max = Pe_max;
        this->thr_prec = thr_prec;
    };
    
    //! Get maximum stable probability mass for degree 2 VN edges
    virtual double get_lam2stable(double sig) = 0;
    
protected:
    //! Ensemble to evolve
    LDPC_Ensemble ens;
    //! Maximum number of iterations
    int maxiter_de;
    //! Maximum number of bisection iterations
    int maxiter_bisec;
    //! Precision of bisection search
    double thr_prec;
    //! Convergence probability of DE
    double Pe_max;
    //! Minimum of bisec threshold search window
    double thr_min;
    //! Maximum of bisec threshold search window
    double thr_max;
    //! Wether the geometric or arithmetic mean should be used for the bisection search
    int mean_mode;
    //! How many de iterations not decreasing the error probability are tolerated before density evolution terminates
    int max_ni_de_iters;
    
    
    //! Wether the object has been initialized
    bool init_flag= false;
};
    
    
/*
 \brief LDPC Density Evolution for purely Lookup-Table (LUT) based decoding
 */
class LDPC_DE_LUT : public LDPC_DE{
    
public:
    //! Default Constructor
    LDPC_DE_LUT(){}
    
    LDPC_DE_LUT(   const LDPC_Ensemble ens_,
                   int Nq_Cha_,
                   ivec Nq_Msg_vec_,
                   int maxiter_de_,
                   Array<Array<LUT_Tree>> var_tree_templates_,
                   Array<Array<LUT_Tree>> chk_tree_templates_,
                   const bvec& reuse_vec_ = bvec(0),
                   double thr_prec_ = 1e-6,
                   double Pe_max_ = 1e-9,
                   int mean_mode_ = ARI,
                   int maxiter_bisec_ = 30,
                   double LLR_max = 25,
                   int Nq_fine_ = 5000,
                   const std::string& irregular_design_strategy_ = "joint_root");
    
    //! Returns the number of iterations if \c thr is feasible for \c ens, otherwise returns -1
    int evolve(double thr, bool var_trace, bool chk_trace, mat& P, vec& p, bool save_luts, Array<Array<LUT_Tree>>& var_trees, Array<Array<LUT_Tree>>& chk_trees);
    
    virtual int evolve(double thr, bool var, bool chk, mat& P, vec& p);
    
    virtual int evolve(double thr);
    
    //! Evolve densities and adaptively set reuse if possible
    bvec evolve_adaptive_reuse(double thr, vec& Pe_trace, double rel_increase_max, double rel_decrease_min, int reuse_max); 
    //! Set the ensemble
    virtual void set_ensemble(const LDPC_Ensemble& ens);
    
    //! Compute the quantized channel pmf
    void set_channel_pmf(double sig);

    //! Return the quantized channel pmf
    vec get_channel_pmf() const{return pmf_cha;}; 

    //! Get maximum stable probability mass for degree 2 VN edges
    virtual double get_lam2stable(double sig);  
    
    //! Get the quantizer boundaries
    void get_quant_bound(double sig, vec& qb_Cha, vec& qb_Msg) const;
    //! Set the resolution of
    void set_tree_template_res();
    
    //! Set the reuse vector
    void set_reuse_vec(const bvec& rv){
        it_assert(rv.length() == maxiter_de, "LDPC_DE_LUT::set_reuse_vec(): Inconsisten length!");
        reuse_vec = rv;
    };
    //! Get the reuse vector 
    bvec get_reuse_vec() const {return reuse_vec;}; 
        /*! \brief Returns an array variable lut tree quantizers designed using density evolution
     
     */
    void get_lut_trees(Array<Array<LUT_Tree>>& var_trees, Array<Array<LUT_Tree>>& chk_trees, double sig);
    
private:
    void chk_update_irr(int iter, Array<LUT_Tree>& prev_trees_chk){
        vec a;
        double b;
        chk_update_irr( a, b, iter, prev_trees_chk);
    }
    void var_update_irr(int iter, Array<LUT_Tree>& prev_trees_var){
        vec a;
        double b;
        var_update_irr( a, b, iter, prev_trees_var);
    }
    void chk_update_irr(vec& P_row, double& Pe, int iter, Array<LUT_Tree>& prev_trees_chk);
    void var_update_irr(vec& P_row, double& Pe, int iter, Array<LUT_Tree>& prev_trees_var);
    
    
    
    
private:
    
    enum{
        INDIVIDUAL,
        JOINT_LEVEL,
        JOINT_ROOT
    };
    
    // Stores a tree structure for each iteration (rows) and variable node degree (columns)
    Array<Array<LUT_Tree>> var_tree_templates;
    Array<Array<LUT_Tree>> chk_tree_templates;
    
    
    vec pmf_cha;
    vec pmf_var2chk;
    vec pmf_chk2var;
    
    bool min_lut;
    bvec reuse_vec;
    ivec Nq_Msg_vec;
    int Nq_Cha;
    //! Number of quantization intervals for fine prequantization used fot getting the channel pmf
    int Nq_fine;
    //! LLR range considered for channel quantization
    double LLR_max;
    //! How the LUTs should be designed in case of irregular ensembles
    int irregular_design_strategy;

    
};

    /*!
     \brief
     Density Evolution for Belief Propagation decoding
     */
class LDPC_DE_BP : public LDPC_DE
{
public:
    LDPC_DE_BP(LDPC_Ensemble ens_,  int Nb_ = 8, double Lmax_ = 25);
    
    virtual int evolve(double thr);
    virtual int evolve(double thr, bool var, bool chk, mat& P, vec& p);
   
    //! Set the ensemble
    virtual void set_ensemble(const LDPC_Ensemble& ens);
   
    //! Get maximum stable probability mass for degree 2 VN edges
    virtual double get_lam2stable(double sig);  
    
private:
    //! Wether the object has been initialized
    
    bool setup_flag = false;
    
    
    
    /*!
     This functions performs a boxplus convolution (c.f. REF) of the input pmf specified by its
     positive and negative parts \c pmf_in_p and \c pmf_in_m and the output pmf, specified
     in a similar manner by \c pmf_out_p and \c pmf_out_m. The result of the convolution
     is again stored in \c pmf_out_p and \c pmf_out_m.
     */
    void chk_update_convolve(const vec& pmf_in_p, const vec& pmf_in_m, vec& pmf_out_p, vec& pmf_out_m);
    
    void chk_update_irr(bool trace, vec& P_row, double& Pe);
    void chk_update_irr(){vec a; double b; chk_update_irr(false,a,b);}
    void var_update_irr(bool trace, vec& P_row, double& Pe);
    void var_update_irr(){vec a; double b; var_update_irr(0,a,b);}
    
    
    
    imat gen_Q_table();
    void set_tq_tables();
    
    
    
    
    
    /*!
     This functions performs a convolution (c.f. REF) of the input pmf \c pmf_in and the output pmf \c pmf_out where the result of the convolution
     is again stored in \c pmf_out.
     */
    void var_update_convolve(const vec& pmf_in, vec& pmf_out);
    
    
    void setup_tables();
    
    /*! \brief Extract the positive part from \c pmf
     
     \c pmf is required to have length 2*N+2, where the first N points are negative, then there is the mass at zero, the next N points are positive and
     the last point is the probability mass at infinity. The output vector has length N+2, as it contains the sum of the onesided parts as well as the
     masses at zero and infinity
     */
    vec pmf_plus(const vec& pmf);
    
    /*! \brief Extract the negative part from \c pmf
     
     \c pmf is required to have length 2*N+2, where the first N points are negative, then there is the mass at zero, the next N points are positive and
     the last point is the probability mass at infinity. The output vector has length N+2, as it contains the sum of the onesided parts as well as  as well as the
     masses at zero (which is 0)  and infinity
     */
    vec pmf_minus(const vec& pmf);
    
    //! Merge positive and negative parts of pmf
    vec pmf_orig(const vec& pmf_p, const vec& pmf_m);
    
    /*!
     Assumes an input vector of length 2*N+2. The output vector is zero padded to length \c Nfft and is
     prepared such that the element at 0 is the first emelment and the element at -delta is the last element
     */
    vec fft_preprocess(const vec& x) const;
    vec ifft_postprocess(const vec& x) const;
    // eliminate negative elements from pmf and scale it to 1
    void numeric_postprocessing(vec& x) const;
    

    vec pmf_LLR;
    vec pmf_var2chk;
    vec pmf_chk2var;

    
    
    //! Maximum LLR magnitude
    double Lmax;
    //! Number of bits
    int Nb;
    //! Total number of pmf points is 2*N+1
    int N;
    //! Difference between x values of two probability mass points
    double delta;
    //! Resolution of FFT
    int Nfft;
    //! Support values of the pmf
    vec support;
    //! Support indices of the pmf (linear)
    ivec idx_support;
    //! Support indices of the pmf (symmetric around 0)
    ivec idx_support_sym;
    //! Transform indices before fft, such that zero element is at first position
    ivec  idx_fft, idx_ifft;
    //! This vector is needed for the Variable node convolution and is precomputed and stored for reasons of efficiency
    vec var_conv_weight;
    //! Tables for check node update
    imat tq, tq2;
    //! Width of check node update tables
    int K;
};



    

   
    

    
vec chk_update_minsum(const vec& p_in, int dc);


    
inline vec pmf_plus(const vec& pmf);
inline vec pmf_minus(const vec& pmf);
inline vec pmf_join(const vec& pmf_p, const vec& pmf_m);

    


    
    

    
void get_lut_tree_templates(const std::string& tree_method, const LDPC_Ensemble& ens, ivec Nq_Msg, int Nq_Cha, bool minLUT, Array<Array<LUT_Tree> >& var_luts, Array<Array<LUT_Tree> >& chk_luts );
    



/*!
 \brief Design trees by taking into account the edge distribution of an irregular LDPC Code
 
 This function takes as input a degree distribution and an array of lut trees.
 The function then aggregates the joint distribution of all lut tree nodes of the same level and performs a joint
 quantizer design for all involved degrees and returns the resulting output pmf. Furthermore, the LUT
 trees are passed as references and are subject to the update of the function.
 */
vec joint_level_irr_lut_design(const vec& degree_dist, const ivec& degrees, Array<LUT_Tree>& lut_trees, vec& P_row, double& Pe);
    
//! Update Trees individually but perform joint design for root luts
vec joint_root_irr_lut_design(const vec& degree_dist, const ivec& degrees, Array<LUT_Tree>& lut_trees, vec& P_row, double& Pe);
    
vec level_lut_tree_update(Array< std::deque<LUT_Tree_Node*> >& tree_nodes,  const vec& degree_dist, LUT_Tree::tree_type_t t);

//! Calculate the maximum stable VN edge degree 2 for a given check node degree, assuming BP decoding on a quantized channel
double get_lam2stable_qbp(double sig, vec rho, int Nq_Cha=5000, double LLR_max=25, int Nq_fine=5000);
//! Iteratively calculate the maximum stable VN edge degree 2 for a given check node degree, assuming BP decoding on a quantized channel
double get_lam2stable_qbp_iterative(double sig, vec rho, int Nq_Cha=pow2i(4), double LLR_max=25, int Nq_fine=pow2i(11));
//! Calculate the maximum stable VN edge degree 2 for a given check node degree, assuming BP decodung on a continous-output channel
double get_lam2stable_cbp(double sig, vec rho);
//! Calculate the maximum stable VN edge degree 2 for a given check node degree, assuming LUT decoding
double get_lam2stable_lut(double sig, vec rho, int Nq_Cha, int Nq_Msg, double LLR_max=25, int Nq_fine=5000);
    


}
#endif
