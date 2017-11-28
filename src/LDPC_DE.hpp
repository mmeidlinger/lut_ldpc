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

// LLRs are considered identical within this tolerance
#define UNIQE_LLR_DELTA 0.0

#ifndef LDPC_DE_hpp
#define LDPC_DE_hpp

#include <deque>
#include <list>
#include <iostream>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iomanip>

//! Git version of program
extern const char *gitversion;

namespace itpp{

// Forward Declatations
class LUT_Tree;
class LUT_Tree_Node;

//! Store degree distribution is sparse form
class LDPC_Ensemble{
private:
    /*! 
     * Check consistency of ensemble.
     */
    inline void check_consistency();
private:
    vec rho;    //< check node degree distrivution (edge perspective) in sparse form
    vec lam;    //< variable node degree distrivution (edge perspective) in sparse form
    ivec degree_rho; //< nonzero check node degrees
    ivec degree_lam; //< nonzero variable node degrees
    int dv_act; //< number of nonzero variable node degrees
    int dc_act; //< number of nonzero check node degrees
    bool init_flag; //< true if the ensemble is defined completely and consistent 
    /*!
     * \brief Accepted deviation to probability mass one without an error
     *
     * When setting any degree distributions, an input check is performed to see wether
     * the degree distributions sum to one, where deviations up to \c pmass_tolerance
     * don't cause an error. Note that the input will still be normalized before assignemnt.
     */
    static const double pmass_tolerance;    
public:
    LDPC_Ensemble();
    //! LDPC ensemble based on nonsparse degree distribution vectors. First element = degree 1
    LDPC_Ensemble(const vec& l, const vec& r);
    //! LDPC ensemble based on sparse degree distribution vectors
    LDPC_Ensemble(const ivec& dl, const vec& l, const ivec& dr, const vec& r);
    //! Read LDPC ensemble from file
    LDPC_Ensemble(const std::string& filename);
    
    //! Read LDPC ensemble from .ens file
    void read(const std::string& filename);
    
    //! Write LDPC ensemble to .ens file
    void write(const std::string& filename) const;
    
    //! Export LDPC ensemble to .deg file
    void export_deg(const std::string& filename) const;
    
    //! Get the rate of the ensemble
    double get_rate() const;
    
    //! Get variable node degree distribution (edge perspective) in sparse form (non zeros only)
    vec sget_lam()const;
    //! Get variable node degree distribution (node perspective) in sparse form (non zeros only)
    vec sget_Lam()const;
    //! Get variable node degree distribution and degrees (edge perspective) in sparse form (non zeros only). Returns the number of active degrees
    int sget_lam(vec& l, ivec& dl)const;
    
    //! Get check node degree distribution (edge perspective) in sparse form (non zeros only)
    vec sget_rho()const;
    //! Get check node degree distribution (node perspective) in sparse form (non zeros only)
    vec sget_Rho()const;
    //! Get check node degree distribution and degrees (edge perspective) in sparse form (non zeros only). Returns the number of active degrees
    int sget_rho(vec& r, ivec& dr)const;
    
    //! Get nonzero variable node degrees
    ivec sget_degree_rho()const;
    //! Get nonzero check node degrees
    ivec sget_degree_lam()const;
    
    //! Get number of nonzero variable node degrees
    int get_dv_act() const;
    //! Get number of nonzero check node degrees
    int get_dc_act() const;
    
    //! Return check node degree distribution (edge perspective) in non-sparse form. First index = degree 1
    vec get_chk_degree_dist() const;
    //! Return variable node degree distribution (edge perspective) in non-sparse form. First index = degree 1
    vec get_var_degree_dist() const;

    //! Set check node degree distribution (edge perspective) in non-sparse form. First index = degree 1
    void set_chk_degree_dist(const vec& r);
    //! Set variable node degree distribution (edge perspective) in non-sparse form. First index = degree 1
    void set_var_degree_dist(const vec& l);
    
    //! Set check node edge distribution in sparse form
    void sset_rho(vec r);
    //! Set variable node edge distribution in sparse form
    void sset_lam(vec l);
    
    //! Get probability mass of VN edge degree \c d 
    double get_lam_of_degree(int d) const;
    //! Print properties of the LDPC ensemble
    friend std::ostream& operator<<(std::ostream &os, const LDPC_Ensemble &ens);
    
};

//! Container class for a tree consisting of LUT_Tree_Nodes
class LUT_Tree{
    public:
        enum tree_type_e{
            VARTREE,
            CHKTREE,
            DECTREE,
            num_tree_types
        }typedef tree_type_t;
        
        
        //! Default constructor, empty tree
        LUT_Tree() : root(nullptr), num_leaves(0) {}
        //! Construct a LUT tree from a string
        LUT_Tree(const std::string& tree_string, tree_type_t t);
        //! Automatically construct a LUT tree of type \c t with \c l leave nodes using method \c m
        LUT_Tree(int l , tree_type_t t , const std::string& m = "auto_bin_balanced");
        //! Copy Constructor
        LUT_Tree(const LUT_Tree& other);
        //! Move Constructor
        LUT_Tree(LUT_Tree && other);
        //! Copy/Move Assignment
        LUT_Tree& operator=(LUT_Tree rhs);
        //! Destructor
        ~LUT_Tree();
        
        
        friend std::ostream& operator<<(std::ostream &os, const LUT_Tree &t);
        friend std::istream& operator>>(std::istream &is, LUT_Tree &t);
        
        //=========== Wrapper functions
        
        //! Returns the cumulated distance of leaf nodes to the tree root
        int get_metric() const;
        //! Set the distribution of the leave nodes
        void set_leaves(const vec& p_Msg, const vec& p_Cha);
        //! Set the resolution of non-root, non-channel nodes to \c Nq_in, the resolution of the root node to \c Nq_out and the resolution of channel LLRs to \c Nq_Cha
        void set_resolution(int Nq_in, int Nq_out, int Nq_cha = 0);
        //! Write tikz code representing the tree to \c outstream
        void tikz_draw_tree(std::ostream& outstream) const;
        //! Write tikz code representing the tree to the file \c filename
        void tikz_draw_tree(const std::string& filename) const;
        //! Get the number of leaf nodes
        int get_num_leaves() const;
        
        //! Output tree structure string
        std::string gen_template_string();
        
        
        //! Get Tree type
        tree_type_t get_type() const {return this->type;}
        /*!
         \brief Updates the non-leaf nodes of the tree by propagating the distribution of the leave nodes to the root via DE
         
         This function is used to update the pmfs of the non-leaf nodes by propagating the distributions of the leave nodes to the root via Density Evolution.
         If \c reuse is true, the quantizers already present within the tree are used for the updates.
         If \c reuse is false, the quantizers are updated as well. The new quangtizers are designed to maximize the local information flow through the tree.
         
         @param[out]    p_out       message pmf at the output of the root node after the update
         @param[out]    mi          mutual information between coded bit and messages at the output of the root node after the update
         @param[in]     reuse       wether new quantizers should be designed or the old ones should be reused
         @param[in]     fp          a function po
         */
        vec update(bool reuse = false);
        
        //! Set the size of the pmfs of all nodes to zero
        void reset_pmfs();
        
        
        //! Variable node message update
        ivec var_msg_update(std::deque<int>& msg_que_all, int llr);
        //! Check node message update
        ivec chk_msg_update(std::deque<int>& msg_que_all);
        //! Decision node llr calculation
        int dec_update(std::deque<int>& msg_que_all, int llr);
        
        
        //! Append pointers to all nodes of the requested level to \c nodes
        void get_level_nodes(int level, std::deque<LUT_Tree_Node*>& nodes);
        
        //! Append pointers to all nodes of the requested level to \c nodes
        std::deque<LUT_Tree_Node*> get_level_nodes(int level);
        
        //! Returns the number of levels of the tree
        int get_height() const;
        
        
    private:
        //! Swaps the objects \c a and \c b
        void swap(LUT_Tree& a, LUT_Tree& b);
        //! Local update of a varaible tree node
        static void var_update(vec& p_out, ivec& Q_out, const Array<vec>& p_in, int Nq, bool reuse);
        //! Locak update if a check tree node
        static void chk_update(vec& p_out, ivec& Q_out, const Array<vec>& p_in, int Nq, bool reuse);
        
        
        
        
    private:
        //! Pointer to root node
        LUT_Tree_Node* root;
        //! Number of tree leaves
        int num_leaves;
        //! Define wether the tree represents a variable or check node update
        tree_type_t type;
    };
    
    //! Print Information about the LDPC ensemble
    std::ostream& operator<<(std::ostream &os, const LDPC_Ensemble &ens);
    
    //Tree I/O
    std::ostream& operator<<(std::ostream &is, const LUT_Tree &t);
    std::istream& operator>>(std::istream &os, LUT_Tree &t);
    
    std::istream& operator>>(std::istream &is, Array<Array<LUT_Tree>> &t);
    std::ostream& operator<<(std::ostream &os, const Array<Array<LUT_Tree>> &t);

//! Super class of all other node nypes
class LUT_Tree_Node{

friend class LUT_Tree;
private:
    enum node_type_e{
        IM,
        ROOT,
        MSG,
        CHA,
        num_node_types
    }typedef node_type_t;
    
    friend std::ostream& operator<<(std::ostream &os, const LUT_Tree &t);
    friend std::istream& operator>>(std::istream &is, LUT_Tree &tree);
    friend vec joint_level_irr_lut_design(const vec& degree_dist,const ivec& degrees, Array<LUT_Tree>& lut_trees, vec& P_row, double& Pe);
    friend vec level_lut_tree_update(Array< std::deque<LUT_Tree_Node*> >& tree_nodes,  const vec& degree_dist, LUT_Tree::tree_type_t t);
    
    //! children of the node
    std::deque<LUT_Tree_Node*> children;
    //! quantizer map
    ivec Q;
    //! pmf
    vec p;
    //! Type of node
    node_type_t type;
    //! number of elements of p and range of Q
    int K;
    
    //! Calls add_child_back()
    void add_child(LUT_Tree_Node* child);
    //! Calls add_child_back()
    void add_child(node_type_t t);
    
    //! Add a child at the end of the children deque
    void add_child_back(LUT_Tree_Node* child);
    //! Create a new node of type \t and add it at the end of the children deque of the calling node
    void add_child_back(node_type_t t);
    
    //! Add a child at the beginning of the children deque
    void add_child_front(LUT_Tree_Node* child);
    //! Create a new node of type \t and add it at the beginning of the children deque of the calling node
    void add_child_front(node_type_t t);
    
    LUT_Tree_Node(node_type_t t);
    
    
    //! Constructor based on input stream
    LUT_Tree_Node(std::istream& is);
    //! Return a string representing the data in the node
    std::string node_to_string();
    
    LUT_Tree_Node* deep_copy();
    int get_metric(int l=0);
    
    //! Set the pmfs of the leave nodes
    void set_leaves(const vec& p_Msg, const vec& p_Cha);
    //! Set the resolution of the tree. ROOT nodes are set to have a resolution of \c Nq_out, CHA nodes \c Nq_cha and all others Nq_in
    void set_resolution(int Nq_in, int Nq_out, int Nq_cha = 0);
    
    //! Returns the number of levels of the tree
    int get_height() const;
    
    vec tree_update( bool reuse,
                     void (*fp)(vec&, ivec&, const Array<vec>&, int, bool));
    /*!
     This function writes  can be saved to a file
     and compiled using e.g., tikz2pdf or with pdflatex if it is
     wrapped appropriately. Most likely, the sibling distance
     needs to be ajusted mannually for the tree nodes not to
     overlap.
     */
    void tikz_draw_tree(std::ostream& outstream);
    void tikz_draw_tree(const std::string& filename);
    void tikz_draw_recursive(std::ostream& outstream, int level=0);
    /*!
     \brief Builds a dynamic LUT tree based on the input string \c s
     */
    static LUT_Tree_Node* parse(std::istream& instream);
    /*!
     \brief This function behaves inverse to parse: It generates a tree template string
     */
    void gen_template_string(std::string& ss);
    
    /*!
     This function returns a tree of \c num_leaves leaves.
     If var == true,
     the tree is created in the following way:  Create binary with dv-1 leaves that is as balanced as possible
     Subsequently attach this tree and a new leave for the channel LLR to a newly created root node 
     If var == false, the attaching of a channel LLR node is skipped.
     */
    static LUT_Tree_Node* gen_bin_balanced_tree(int num_leaves, bool var, node_type_t leaf_type = MSG);
    
    /*!
     This function returns a tree of \c num_leaves leaves.
     The tree is binary and has maximum height, thus resembling the approach of 
        Lewandowsky et al.: Trellis based node operations for LDPC decoders from the Information Bottleneck method
        Romero et al.: Decoding LDPC codes with mutual information-maximizing lookup tables
        Kurkoski et al.: Noise Thresholds for Discrete LDPC Decoding Mappings
     */
    static LUT_Tree_Node* gen_bin_high_tree(int num_leaves, bool var, node_type_t leaf_type = MSG);
    
    /*!
        Generate a tree that only consists of a root node with all messages/ LLRs directly attached to it
    */
    static LUT_Tree_Node* gen_root_only_tree(int num_leaves, bool var, node_type_t leaf_type = MSG);
    
    //! Delete all descendant nodes
    void delete_tree();
    
    //! Traverse nodes to free the memory allocated for pmfs
    void reset_pmfs();
    
    //! Return the number of leaf nodes
    int get_num_leaves();
    
    //! Return the number of LUT nodes
    int get_num_luts();
    
    //! Calculate LUT tree output for input messages \c msgs incident to the leaf nodes
    int var_msg_update(std::deque<int>& msgs);
    
    //! Calculate LUT tree output for input messages \c msgs incident to the leaf nodes
    int chk_msg_update(std::deque<int>& msgs);

    /*!
    \brief Serialize the tree recursively
    */
    void serialize_tree(std::ostream& os);
    /*!
    \brief Deserialize the tree recursively
    */
    static LUT_Tree_Node* deserialize_tree(std::istream& is);
    
    //! get the product pmf of the nodes children
    vec get_input_product_pmf(LUT_Tree::tree_type_t t) const;
    
    //! return pointers to all nodes (\c level_nodes) of a certain level \c req_level
    void get_level_nodes(int req_level, int cur_level, std::deque<LUT_Tree_Node*>& level_nodes);
};

    
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



    
/*!
 * Compute the mutual information optimal quantizer for the symmetric pmf p_in. Returns the resulting Mutual information
 *
 * @param[out] p_out    Conditional length Nq output pmf
 * @param[out] Q_out    Quantizer designed to maximize the mutual informatio \c mi. Implementation wise, this is a
 *                      length M vectors with integer elements in 0,1,...,\c Nq-1
 * @param[in]  p_in     Conditional input pmf of length M. It is assumed, that this is a symmetric conditional pmf, where the conditional random variable X is uniform and binary:
 *                      i.e., p(y|x) = p(-y|-x), where p_in[M/2+i-1] = p(i|1), i=-M/2,...,M/2
 * @param[in]  Nq       Number of quantizer outputs
 * @param[in]  Q_old    If the quanzizer should be reused rather than designed, this is indicated by passing a vector \c Q_old with length > 0
 */
double quant_mi_sym(vec& p_out, ivec& Q_out, const vec& p_in, int Nq, bool sorted = false);
   
    
/*!
 * \brief
 * Returns a symmetric pmf with unique LLRs and the corresponding indices.
 *
 *
 * @returns	 Returns a symmetric pmf with unique LLRs and the corresponding indices
 * @param[in]   p_in        An arbitrary, unsorted, conditional pmf.
 * @param[out]  idx_in      sorting index of input pmf according to LLR
 * @param[out]  idx_sorted  index mapping input to unique pmf
 */
vec sym_llr_sort_unique(const vec& p_in, ivec& idx_in, ivec& idx_sorted, double llr_delta = UNIQE_LLR_DELTA);
    
vec chk_update_minsum(const vec& p_in, int dc);

//! Calculate the mutual information between X and Y, where p(y|x)=p(-y|-x) is given by p_in and X is binary and uniform
double get_mi_bcpmf_sym(const vec& p);
    
inline vec pmf_plus(const vec& pmf);
inline vec pmf_minus(const vec& pmf);
inline vec pmf_join(const vec& pmf_p, const vec& pmf_m);

    

int quant_nonlin(double x, const vec& boundaries);
int quant_lin(double x, double delta, int N);
ivec quant_nonlin(const vec& x, const vec& boundaries);
/*! \brief Get quantized Gaussian pmfs with N quantization intervals (2 Overload and
 N -2) inner regions. For N odd, this is quantization with 0, with N even there is no zero.
 */
vec get_gaussian_pmf(double mu, double sig, int N, double delta);

double rate_to_shannon_thr(double R);
double shannon_thr_to_rate(double sig);
    

inline double x_log2_y(double x, double y);
    
template <class Num_T> Vec<Num_T> fliplr(const Vec<Num_T>& x);
template <class Num_T> Vec<Num_T> unique(const Vec<Num_T>& x);
template <class Num_T> Vec<Num_T> kron(const Vec<Num_T>& x, const Vec<Num_T>& y);
    
LDPC_Ensemble get_empirical_ensemble(const LDPC_Parity& H);
    
template <class Num_T> ivec sort_index_sym(const Vec<Num_T>& x);
template<class Num_T> ivec get_complement_idx(const Vec<Num_T>& a, const Vec<Num_T>& b,  Num_T c);
    
int signed_to_unsigned_idx(int idx, const ivec& inres);
    
void get_lut_tree_templates(const std::string& tree_method, const LDPC_Ensemble& ens, ivec Nq_Msg, int Nq_Cha, bool minLUT, Array<Array<LUT_Tree> >& var_luts, Array<Array<LUT_Tree> >& chk_luts );
    

//! Get the product distribution for a variable node with input probabilities p_in
vec get_var_product_pmf(const Array<vec>& p_in);

//! Get the product distribution for a check node with input probabilities p_in
vec get_chk_product_pmf(const Array<vec>& p_in);

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
    
//! Convert SNR (=Eb/N0 in dB) to corresponding AWGN channel noise standard deviation
inline double snr2sig(double rate, double snr);
//! Convert AWGN channel noise standard deviation to SNR (=Eb/N0 in dB) 
inline double sig2snr(double rate, double sig);

//! Convert AWGN channel noise standard deviation to SNR (=Eb/N0 in dB) 
vec sig2snr(double rate, const vec& sig);

//! Convert SNR (=Eb/N0 in dB) to corresponding AWGN channel noise standard deviation
vec snr2sig(double rate, const vec& snr);

}
#endif
