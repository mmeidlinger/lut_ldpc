/*!
 * \file LUT_Tree.hpp
 * \brief Implementation of tree structured Lookup Table (LUT) node updates for discrete LDPC message passing decoding
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

#ifndef LUT_Tree_hpp
#define LUT_Tree_hpp

#include "LDPC_Ensemble.hpp"
#include "common.hpp"
#include <stdio.h>
#include <itpp/itbase.h>
#include <deque>
#include <list>

using namespace itpp;

namespace lut_ldpc{
    // Forward Declatations
    class LUT_Tree;
    class LUT_Tree_Node;
    
    //! Container class for a tree consisting of LUT_Tree_Nodes
    class LUT_Tree{
    public:
        //! Type of tree, can be either variable node, check node or decision node
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
        friend std::istream& operator>>(std::istream &is, Array<Array<LUT_Tree>> &t);
        friend std::ostream& operator<<(std::ostream &os, const Array<Array<LUT_Tree>> &t);

        
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
    /*!
        \brief Serialize tree and write it to output stream
     */
    std::ostream& operator<<(std::ostream &is, const LUT_Tree &t);

    /*!
     \brief Construct tree by reading a serilized version from input stream
     */
    std::istream& operator>>(std::istream &os, LUT_Tree &t);
    
    /*!
     \brief Serialize array of trees and write them to output stream
     
     The first line contains the number of iterations (first array dimension), referred to as I. Subsequently, for each iteration,
     the stream contains one tree for each of the D degrees (second array dimension). The LUT trees themselves are serialized
     via LUT_Tree_Node::serialize_tree()
     */
    std::istream& operator>>(std::istream &is, Array<Array<LUT_Tree>> &t);
    
    /*!
     \brief Construct array of trees by reading a serilized version from input stream
     */
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
        
        //! Deep copy the tree with all its nodes
        LUT_Tree_Node* deep_copy();
        
        //! Get the cumulative distance from leave nodes to the root node
        int get_metric(int l=0);
        
        //! Set the pmfs of the leave nodes
        void set_leaves(const vec& p_Msg, const vec& p_Cha);
        //! Set the resolution of the tree. ROOT nodes are set to have a resolution of \c Nq_out, CHA nodes \c Nq_cha and all others Nq_in
        void set_resolution(int Nq_in, int Nq_out, int Nq_cha = 0);
        
        //! Returns the number of levels of the tree
        int get_height() const;
        
        //! Updates the LUTs recursively, calling the function pointer \c fp for each tree node. Depending on the type of tree, \c fp would be set to point to LUT_Tree::var_update() or LUT_Tree::chk_update()
        vec tree_update( bool reuse,
                        void (*fp)(vec&, ivec&, const Array<vec>&, int, bool));
        /*!
         This function serves to visualize the tree structure. generates TikZ LaTeX code which can be
         compiled with tikz2pdf or with pdflatex if it is
         wrapped appropriately. Most likely, the sibling distance
         needs to be ajusted mannually for the tree nodes not to
         overlap.
         */
        void tikz_draw_tree(std::ostream& outstream);
        /*!
            Same as tikz_draw_tree(std::ostream& outstream), except the output is written to file
         */
        void tikz_draw_tree(const std::string& filename);
        
        /*!
         Actual functionality behind tikz_draw_tree(std::ostream& outstream), drawing tree nodes recursively
         */
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
         \brief Serialize the tree recursively, using in order depth-first traversation
         */
        void serialize_tree(std::ostream& os);
        /*!
         \brief Deserialize the tree recursively, using in order depth-first traversation
         */
        static LUT_Tree_Node* deserialize_tree(std::istream& is);
        
        //! get the product pmf of the nodes children
        vec get_input_product_pmf(LUT_Tree::tree_type_t t) const;
        
        //! return pointers to all nodes (\c level_nodes) of a certain level \c req_level
        void get_level_nodes(int req_level, int cur_level, std::deque<LUT_Tree_Node*>& level_nodes);
    };
}
#endif /* LUT_Tree_hpp */
