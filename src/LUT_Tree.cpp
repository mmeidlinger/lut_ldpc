/*!
 * \file LUT_Tree.cpp
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

#include "LUT_Tree.hpp"
using namespace lut_ldpc;
using namespace itpp;
using namespace std;

// ----------------------------------------------------------------------
// LUT_Tree_Node
// ----------------------------------------------------------------------

LUT_Tree_Node::LUT_Tree_Node(node_type_t t){
    type = t;
    K = 0;
    ivec Q(0);
    vec p(0);
}

void LUT_Tree_Node::add_child(node_type_t t){
    this->add_child_back(t);
    return;
}

void LUT_Tree_Node::add_child(LUT_Tree_Node* child){
    this->add_child_back(child);
    return;
}

void LUT_Tree_Node::add_child_front(node_type_t t){
    LUT_Tree_Node* new_child = new LUT_Tree_Node(t);
    children.push_front(new_child);
    return;
}

void LUT_Tree_Node::add_child_front(LUT_Tree_Node* child){
    if( child != 0)
        children.push_front(child);
    return;
}

void LUT_Tree_Node::add_child_back(node_type_t t){
    LUT_Tree_Node* new_child = new LUT_Tree_Node(t);
    children.push_back(new_child);
    return;
}

void LUT_Tree_Node::add_child_back(LUT_Tree_Node* child){
    if( child != 0)
        children.push_back(child);
    return;
}

LUT_Tree_Node* LUT_Tree_Node::deep_copy(){
    
    LUT_Tree_Node* new_root =  new LUT_Tree_Node(this->type);
    new_root->p = this->p;
    new_root->Q = this->Q;
    new_root->K = this->K;
    
    for (unsigned int ii=0; ii< this->children.size(); ii++){
        new_root->children.push_back( this->children[ii]->deep_copy());
    }
    return new_root;
}

void LUT_Tree_Node::set_leaves(const vec& p_Msg, const vec& p_Cha){
    if(type == MSG)
        p = p_Msg;
    else if(type == CHA)
        p = p_Cha;
    else{
        for (unsigned int ii=0; ii<children.size(); ii++){
            children[ii]->set_leaves(p_Msg,p_Cha);
        }
    }
    return;
}

int LUT_Tree_Node::get_height() const{
    int h = 0;
    for(unsigned int ii=0; ii<children.size(); ii++){
        int h_tmp = children[ii]->get_height();
        if( h_tmp >= h)  h = h_tmp+1;
    }
    return h;
}

vec LUT_Tree_Node::tree_update(bool reuse,
                               void (*fp)(vec&, ivec&, const Array<vec>&, int, bool) )
{
    if( type == MSG || type == CHA){
        return p;;
    }
    else{
        Array<vec> p_Msg((int)children.size());
        for(unsigned int ii=0; ii<children.size(); ii++){
            p_Msg(ii) = children[ii]->tree_update(reuse, fp);
        }
        // At this point we are at an im or root node and and have examined and collected
        // all children pmfs in p_Msg
        (*fp)(p, Q, p_Msg, K, reuse);
        return p;
    }
}

int LUT_Tree_Node::get_metric(int l){
    if(this->type == MSG || this->type == CHA)
        return l;
    else{
        l++;
        for(unsigned int ii=0; ii<children.size(); ii++)   l += children[ii]->get_metric();
    }
    return l;
}

void LUT_Tree_Node::gen_template_string(std::string& ss){
    
    switch (this->type) {
        case ROOT:
            ss =  ss + "r";
            break;
        case IM:
            ss =  ss + "i";
            break;
        case MSG:
            ss =  ss + "m";
            break;
        case CHA:
            ss =  ss + "c";
            break;
        default:
            break;
    }
    
    for(unsigned int ii=0; ii<children.size(); ii++)
        this->children[ii]->gen_template_string(ss);
    
    ss = ss + "/";
}

LUT_Tree_Node* LUT_Tree_Node::parse(std::istream& instream){
    LUT_Tree_Node* node = 0;
    LUT_Tree_Node* child;
    char c = instream.get();
    switch (c) {
        case EOF:
            return node;
        case '/':
            return node;
        case 'r':
            node = new LUT_Tree_Node(ROOT);
            break;
        case 'i':
            node = new LUT_Tree_Node(IM);
            break;
        case 'm':
            node = new LUT_Tree_Node(MSG);
            break;
        case 'c':
            node = new LUT_Tree_Node(CHA);
            break;
        default:
            it_error("LUT_Tree_Node::parse(): Allowed characters are r (root), i (intermediate node), m (message node), c (channel node) and / (end of children)");
            break;
    }
    while(1){
        child = parse(instream);
        if(child ==0)    break;
        node->add_child(child);
    }
    return node;
}

LUT_Tree_Node* LUT_Tree_Node::gen_bin_balanced_tree(int num_leaves, bool var, node_type_t leaf_type){
    
    it_assert(num_leaves>=2, "LUT_Tree_Node::gen_bin_balanced_tree(): num_leaves must be larger than or equal to 2");
    std::list<LUT_Tree_Node*> nodes;
    for(int ll=0; ll< num_leaves - (int)var; ll++)
        nodes.push_back(new LUT_Tree_Node(leaf_type));
    
    
    LUT_Tree_Node* left;
    LUT_Tree_Node* right;
    LUT_Tree_Node* new_node=0;
    while(true){
        if(nodes.size()==1){
            if(var){
                new_node = new LUT_Tree_Node(ROOT);
                new_node->add_child(nodes.front());
                new_node->add_child(CHA);
                nodes.clear();
                break;
            }
            else{
                new_node = nodes.front();
                nodes.clear();
                new_node->type = ROOT;
                break;
            }
        }
        left = nodes.front();
        nodes.pop_front();
        right = nodes.front();
        nodes.pop_front();
        new_node = new LUT_Tree_Node(IM);
        new_node->add_child(left);
        new_node->add_child(right);
        nodes.push_back(new_node);
    }
    return new_node;
}


LUT_Tree_Node* LUT_Tree_Node::gen_bin_high_tree(int num_leaves, bool var, node_type_t leaf_type){
    
    it_assert(num_leaves>=2, "LUT_Tree_Node::gen_bin_high_tree(): num_leaves must be larger than or equal to 2");
    LUT_Tree_Node* root;
    LUT_Tree_Node* cur_node;
    
    // Add a root node and an intermediate node if necesary
    root = new LUT_Tree_Node(ROOT);
    cur_node = root;
    
    if(var){
        cur_node->add_child(CHA);
    }
    else{
        cur_node->add_child(leaf_type);
    }
    int num_leaves_todo = num_leaves-1;
    
    // Add all but the last leaf nodes
    while(num_leaves_todo>1){
        cur_node->add_child_front(IM);
        cur_node = cur_node->children.front();
        cur_node->add_child(leaf_type);
        num_leaves_todo--;
    }
    
    // Add the last leaf node
    cur_node->add_child(leaf_type);
    
    return root;
}

LUT_Tree_Node* LUT_Tree_Node::gen_root_only_tree(int num_leaves, bool var, node_type_t leaf_type){
    
    it_assert(num_leaves>=2, "LUT_Tree_Node::gen_root_only_tree(): num_leaves must be larger than or equal to 2");
    LUT_Tree_Node* root;
    
    // Add a root node and an intermediate node if necesary
    root = new LUT_Tree_Node(ROOT);
    
    
    for(int ii=0; ii<num_leaves-1; ii++){
        root->add_child_back(leaf_type);
    }
    
    if(var){
        root->add_child_back(LUT_Tree_Node::CHA);
    }
    else{
        root->add_child_back(leaf_type);
    }
    
    
    return root;
}

void LUT_Tree_Node::set_resolution(int Nq_in, int Nq_out, int Nq_cha){
    if (type == ROOT)
        K = Nq_out;
    else if (type == CHA)
        K = Nq_cha;
    else
        K = Nq_in;
    
    for (unsigned int ii=0; ii< children.size(); ii++)
        children[ii]->set_resolution(Nq_in, Nq_out, Nq_cha);
}

void LUT_Tree_Node::tikz_draw_tree(std::ostream& outstream){
    outstream <<  "\\tikzset{\n" <<
    "   leavenode/.style = {align=center, inner sep=2pt, text centered },\n" <<
    "   imnode/.style = {align=center, inner sep=1pt, text centered},\n" ;
    
    int height = this->get_height();
    for(int hh=1; hh<= height; hh++ )
        outstream << "   level " << hh << "/.style={sibling distance=" << 7*pow2i(height-hh) << "mm},\n";
    
    outstream <<  "}\n\n" <<
    "\\def\\imstring{$\\Phi$}\n" <<
    "\\def\\chastring{$L$}\n" <<
    "\\def\\msgstring{$\\mu$}\n\n" <<
    "\\begin{tikzpicture}[<-, >=stealth]";
    
    this->tikz_draw_recursive(outstream);
    
    outstream << "\n\\end{tikzpicture}";
}

void LUT_Tree_Node::tikz_draw_tree(const std::string& filename){
    std::ofstream f;
    f.open(filename);
    this->tikz_draw_tree(f);
    f.close();
}

void LUT_Tree_Node::tikz_draw_recursive(std::ostream& outstream, int level){
    
    // write a new line and indetation according to the current level
    outstream << std::endl;
    for(int ii=0; ii<level; ii++)   outstream << "   ";
    // write nodes according to type
    switch (this->type) {
        case ROOT:
            outstream << "\\node (root)[imnode] {\\imstring}";
            break;
        case MSG:
            outstream << "child{ node [leavenode] {\\msgstring}";
            break;
        case CHA:
            outstream << "child{ node [leavenode] {\\chastring}";
            break;
        case IM:
            outstream << "child{ node[imnode] {\\imstring}";
            break;
        default:
            it_error("LUT_Tree_Node::tikz_draw_recursive(): Node type undefined");
            break;
    }
    // proceed to traverse tree
    for(unsigned int ii=0; ii<this->children.size(); ii++){
        this->children[ii]->tikz_draw_recursive(outstream, level+1);
    }
    // Befor returning, the child{... statements need to be closed
    outstream << std::endl;
    for(int ii=0; ii<level; ii++)   outstream << "   ";
    if(this->type == ROOT)  outstream << ";";
    else    outstream << "}";
    return;
}

void LUT_Tree_Node::delete_tree(){
    for(unsigned int ii=0; ii<this->children.size(); ii++){
        this->children[ii]->delete_tree();
    }
    delete this;
    return;
}

int LUT_Tree_Node::get_num_leaves(){
    if(this->type == MSG || this->type == CHA)      return 1;
    else{
        int nl = 0;
        
        for(unsigned int ii=0; ii<this->children.size(); ii++){
            nl += this->children[ii]->get_num_leaves();
        }
        return nl;
    }
}

int LUT_Tree_Node::get_num_luts(){
    it_error("LUT_Tree_Node::num_luts(): Not implemented yet!");
    return 0;
}

void LUT_Tree_Node::reset_pmfs(){
    for(unsigned int ii=0; ii<this->children.size(); ii++)
        this->children[ii]->reset_pmfs();
    this->p.set_size(0);
}


int LUT_Tree_Node::var_msg_update(std::deque<int>& msgs){
    if(this->type == MSG || this->type == CHA){
        int out = msgs.front();
        msgs.pop_front();
        return out;
    }
    int label = 0; // LUT input
    int base = 1;
    for(unsigned int ii=0; ii<this->children.size(); ii++){
        label += base*this->children[ii]->var_msg_update(msgs);
        base  *= this->children[ii]->K;
    }
    if(label < length(Q))
        return Q(label);
    else
        return K-1 - Q(2*length(Q)-1-label);
}

int LUT_Tree_Node::chk_msg_update(std::deque<int>& msgs){
    if(this->type == MSG ){
        int out = msgs.front();
        msgs.pop_front();
        return out;
    }
    int label = 0; // LUT input
    int base = 1;
    int parity = 0;
    for(unsigned int ii=0; ii<this->children.size(); ii++){
        int label_signed = this->children[ii]->chk_msg_update(msgs);
        int child_res = this->children[ii]->K;
        if( label_signed < child_res/2){
            parity ^= 1;
            label += base* (child_res/2 -1 -label_signed);
        }
        else{
            label += base* (label_signed- child_res/2);
        }
        base  *= child_res/2;
    }
    if(parity == 1)
        return Q(label);
    else
        return K-1 - Q(label);
}


LUT_Tree_Node::LUT_Tree_Node(std::istream& is){
    
    std::string line;
    std::stringstream s;
    int t, inres, outres;
    
    // Get Node Type, Input Resolution and Output Resolution
    getline(is, line);
    s << line;
    s >> t >> inres >> outres;
    it_assert(!s.fail(),
              "LUT_Tree_Node::LUT_Tree_Node(): Error reading node type and resolution");
    it_assert((t >= 0 && t < LUT_Tree_Node::num_node_types) && (inres >= 0) && (outres >= 0),
              "LUT_Tree_Node::LUT_Tree_Node(): Wrong node data");
    s.seekg(0, std::ios::end);
    s.clear();
    
    this->type = (node_type_t) t;
    this->K = outres;
    
    // Parse Quantizer if present
    if(inres > 0){
        getline(is, line);
        s << line;
        ivec map(inres);
        for(int ii=0; ii< inres; ii++){
            s >> map(ii);
            it_assert(!s.fail(),
                      "LUT_Tree_Node::LUT_Tree_Node(): Wrong mapping data (map("
                      << ii << "))");
            it_assert( (map(ii)>=0) && ( map(ii) < K ),
                      "LUT_Tree_Node::LUT_Tree_Node(): Wrong mapping data (map("
                      << ii << "))");
        }
        this->Q = map;
        
    }
}


std::string LUT_Tree_Node::node_to_string(){
    std::stringstream s;
    
    int inres = this->Q.length();
    // First Line: type, inres, outres
    s << static_cast<int>(this->type) << " " << inres << " " << this->K << std::endl ;
    // Second Line: Quanzier Mapping
    if(inres > 0){
        for(int ii=0; ii< inres-1; ii++){
            s << this->Q(ii) << " ";
        }
        s << Q(inres-1) << std::endl;
    }
    return s.str();
}




void LUT_Tree_Node::serialize_tree(std::ostream& os){
    int num_children = (int) this->children.size();
    
    os << num_children << std::endl;
    os << this->node_to_string();
    
    for(int ii=0; ii<num_children; ii++){
        this->children[ii]->serialize_tree(os);
    }
    
}


LUT_Tree_Node* LUT_Tree_Node::deserialize_tree(std::istream& is){
    
    int num_children;
    std::string line;
    std::stringstream s;
    getline(is, line);
    s << line;
    s >> num_children;
    
    LUT_Tree_Node* new_root =  new LUT_Tree_Node(is);
    
    for (int ii=0; ii< num_children; ii++){
        new_root->children.push_back( deserialize_tree(is));
    }
    return new_root;
}

vec LUT_Tree_Node::get_input_product_pmf(LUT_Tree::tree_type_t t) const
{
    
    it_assert(this->type == IM || this->type == ROOT, "LUT_Tree::get_input_product_pmf(): function only supported for IM and ROOT nodes!");
    // fetch child distributions
    int Nchildren =  (int) this->children.size();
    Array<vec> p_in(Nchildren);
    for(int ii=0; ii<Nchildren; ii++){
        p_in(ii) = this->children[ii]->p;
    }
    // build joint distribution
    if( t == LUT_Tree::VARTREE){
        return get_var_product_pmf(p_in);
    }
    else if( t == LUT_Tree::CHKTREE){
        return get_chk_product_pmf(p_in);
    }
    else if( t == LUT_Tree::DECTREE){
        return get_var_product_pmf(p_in);
    }
    else{
        it_error("LUT_Tree::get_input_product_pmf(): Function not implemented for this type of tree");
        return vec(0);
    }
}


void LUT_Tree_Node::get_level_nodes(int req_level, int cur_level, std::deque<LUT_Tree_Node*>& level_nodes){
    if(req_level == cur_level){
        level_nodes.push_back(this);
    }
    else{
        for(unsigned int ii=0; ii<this->children.size(); ii++)
            this->children[ii]->get_level_nodes(req_level, cur_level+1, level_nodes);
    }
}


// ----------------------------------------------------------------------
// LUT_Tree
// ----------------------------------------------------------------------

LUT_Tree::LUT_Tree(const std::string& tree_string, tree_type_t t){
    // check wether the tree string represents a valid variable or check node update
    if( tree_string.find("c") == std::string::npos && t != CHKTREE ){
        it_error("LUT_Tree::LUT_Tree(): Trees other than type CHKTREE need to have a channel leaf specified");
    }
    this->type =t;
    
    // create tree
    std::stringstream tree_stream;
    tree_stream << tree_string;
    
    this->root = LUT_Tree_Node::parse(tree_stream);
    this->num_leaves = root->get_num_leaves();
}

LUT_Tree::LUT_Tree(int l , tree_type_t t  , const std::string& m){
    switch (t) {
        case DECTREE: //same as VARTREE
        case VARTREE:
            if(m == "auto_bin_balanced")
                this->root = LUT_Tree_Node::gen_bin_balanced_tree(l, true);
            else if(m == "auto_bin_high")
                this->root = LUT_Tree_Node::gen_bin_high_tree(l, true);
            else if(m == "root_only")
                this->root = LUT_Tree_Node::gen_root_only_tree(l, true);
            else
                it_error("LUT_Tree::LUT_Tree(): Autogeneration mode " << m << " not supported");
            
            this->type = t;
            this->num_leaves = l;
            break;
            
        case CHKTREE:
            if(m == "auto_bin_balanced")
                this->root = LUT_Tree_Node::gen_bin_balanced_tree(l, false);
            else if(m == "auto_bin_high")
                this->root = LUT_Tree_Node::gen_bin_high_tree(l, false);
            else if(m == "root_only")
                this->root = LUT_Tree_Node::gen_root_only_tree(l, false);
            else
                it_error("LUT_Tree::LUT_Tree(): Autogeneration mode " << m << " not supported");
            
            this->type = t;
            this->num_leaves = l;
            break;
            
            
        default:
            it_error("LUT_Tree::LUT_Tree(): Type must be VARTREE, DECTREE or CHKTREE");
            break;
    }
}

LUT_Tree::LUT_Tree(const LUT_Tree& other){
    if(other.root == nullptr)
        this->root = nullptr;
    else
        this->root = other.root->deep_copy();
    this->num_leaves = other.num_leaves;
    this->type = other.type;
}

LUT_Tree::LUT_Tree(LUT_Tree&& other){
    this->root = other.root;
    other.root = nullptr;
    this->num_leaves = std::move(other.num_leaves);
    this->type = std::move(other.type);
}

LUT_Tree& LUT_Tree::operator=(LUT_Tree rhs){
    LUT_Tree::swap(*this, rhs);
    return *this;
}


LUT_Tree::~LUT_Tree(){
    if(this->root != nullptr)
        this->root->delete_tree();
}


int LUT_Tree::get_metric() const{
    return this->root->get_metric();
}

void LUT_Tree::set_leaves(const vec& p_Msg, const vec& p_Cha){
    this->root->set_leaves(p_Msg, p_Cha);
}

void LUT_Tree::set_resolution(int Nq_in, int Nq_out, int Nq_cha){
    this->root->set_resolution(Nq_in, Nq_out, Nq_cha);
}

void LUT_Tree::tikz_draw_tree(std::ostream& outstream) const {
    this->root->tikz_draw_tree(outstream);
}
void LUT_Tree::tikz_draw_tree(const std::string& filename) const{
    this->root->tikz_draw_tree(filename);
}

int LUT_Tree::get_num_leaves() const{
    return this->num_leaves;
}

vec LUT_Tree::update(bool reuse){
    switch (this->type) {
        case DECTREE: // Same as Var tree, except reuse is prohibited
            // if(reuse){  it_error("Reuse of Decision node tree not allowed"); }
        case VARTREE:
            return this->root->tree_update(reuse, &LUT_Tree::var_update);
            break;
        case CHKTREE:
            return this->root->tree_update(reuse, &LUT_Tree::chk_update);
            break;
        default:
            it_error("LUT_Tree::update(): Tree type must be either VARTREE or CHKTREE");
            return vec(0);
            break;
    }
}


void LUT_Tree::swap(LUT_Tree& a, LUT_Tree& b){
    std::swap(a.type, b.type);
    std::swap(a.root, b.root);
    std::swap(a.num_leaves, b.num_leaves);
    return;
}


void LUT_Tree::var_update(vec& p_out, ivec& Q_out, const Array<vec>& p_in, int Nq, bool reuse){
    
    vec p_Msg_prod = get_var_product_pmf(p_in);
    int M = length(p_Msg_prod);
    
    // Get quantizer, either by reususe or design
    if(reuse){
        p_out = zeros(Nq);
        for(int mm=0; mm<M; mm++){
            if(mm<M/2)
                p_out(Q_out(mm)) += p_Msg_prod(mm);
            else
                p_out(Nq-1-Q_out(M-1-mm)) += p_Msg_prod(mm);
        }
    }
    else{
        // Eliminate entries with mass 0
        bvec idx_nz = (.5*(p_Msg_prod + fliplr(p_Msg_prod)) != 0);
        ivec Q_out_nz;
        (void) quant_mi_sym(p_out, Q_out_nz, p_Msg_prod.get(idx_nz), Nq);
        // For the entries with mass zero we set the outputs symmetrically and assign the least confident llr magnitudes
        Q_out = concat(ones_i(M/2)*(Nq/2-1), ones_i(M/2)*(Nq/2) );
        int mm_nz=0;
        for(int mm=0; mm<M; mm++){
            if(idx_nz(mm)){
                Q_out(mm) = Q_out_nz(mm_nz);
                mm_nz++;
            }
        }
        Q_out = Q_out.left(length(Q_out)/2);
    }
    p_out = p_out/sum(p_out);
    return;
}

void LUT_Tree::chk_update(vec& p_out, ivec& Q_out, const Array<vec>& p_in, int Nq, bool reuse){
    
    vec p_Msg_prod_combined = get_chk_product_pmf(p_in);
    
    // Get quantizer, either by reususe or design
    if(reuse){
        p_out = zeros(Nq);
        int M = length(p_Msg_prod_combined);
        for(int mm=0; mm<M; mm++){
            if(mm<M/2)
                p_out(Q_out(mm)) += p_Msg_prod_combined(mm);
            else
                p_out(Nq-1-Q_out(M-1-mm)) += p_Msg_prod_combined(mm);
        }
        
    }
    else{
        (void) quant_mi_sym(p_out, Q_out, p_Msg_prod_combined, Nq);
        Q_out =  Q_out.left(length(Q_out)/2);
    }
    p_out = p_out/sum(p_out);
    return;
}

void LUT_Tree::reset_pmfs(){
    if(this->root != nullptr)   root->reset_pmfs();
}



ivec LUT_Tree::var_msg_update(std::deque<int>& msg_que_all, int llr){
    it_assert_debug(this->type == VARTREE, "LUT_Tree::var_msg_update(): Only supported for variable node trees");
    it_assert_debug((int) msg_que_all.size() == this->num_leaves, "LUT_Tree::var_msg_update(): Number of inputs must match number of leaves");
    it_assert_debug(this->root != nullptr, "LUT_Tree::var_msg_update(): Tree is empty");
    // todo: check range
    
    ivec msgs_out((int) msg_que_all.size());
    msg_que_all.push_back(llr);
    
    for(int ii=0; ii<length(msgs_out); ii++){
        std::deque<int> msg_que = msg_que_all;
        msg_que.erase(msg_que.begin()+ii);
        msgs_out(ii) = root->var_msg_update(msg_que);
        it_assert_debug(msg_que.size() == 0, "LUT_Tree::var_msg_update(): Not all messages have been processed");
    }
    return msgs_out;
}

ivec LUT_Tree::chk_msg_update(std::deque<int>& msg_que_all){
    it_assert_debug(this->type == CHKTREE, "LUT_Tree::chk_msg_update(): Only supported for check node trees");
    it_assert_debug((int) msg_que_all.size() == this->num_leaves+1, "LUT_Tree::chk_msg_update(): Number of inputs must match number of leaves");
    it_assert_debug(this->root != nullptr, "LUT_Tree::chk_msg_update(): Tree is empty");
    // todo: check range
    
    ivec msgs_out((int) msg_que_all.size());
    
    for(int ii=0; ii<length(msgs_out); ii++){
        std::deque<int> msg_que = msg_que_all;
        msg_que.erase(msg_que.begin()+ii);
        msgs_out(ii) = root->chk_msg_update(msg_que);
        it_assert_debug(msg_que.size() == 0, "LUT_Tree::chk_msg_update(): Not all messages have been processed");
    }
    return msgs_out;
}

int LUT_Tree::dec_update(std::deque<int>& msg_que_all, int llr){
    it_assert_debug(this->type == DECTREE, "LUT_Tree::dec_update(): Only supported for decision node trees");
    it_assert_debug((int) msg_que_all.size() + 1 == this->num_leaves, "LUT_Tree::dec_update(): Number of inputs including llr must match number of leaves");
    it_assert_debug(this->root != nullptr, "LUT_Tree::dec_update(): Tree is empty");
    // todo: check range
    
    msg_que_all.push_back(llr);
    int llr_out = root->var_msg_update(msg_que_all);
    it_assert_debug(msg_que_all.size() == 0, "LUT_Tree::dec_update(): Not all messages have been processed");
    
    return llr_out;
}

std::string LUT_Tree::gen_template_string(){
    std::string ss;
    this->root->gen_template_string(ss);
    return ss;
}


int LUT_Tree::get_height() const
{
    return root->get_height();
}

void LUT_Tree::get_level_nodes(int level, std::deque<LUT_Tree_Node*>& nodes){
    it_assert(root != nullptr, "LUT_Tree::get_level_nodes(): Tree empty!");
    root->get_level_nodes(level, 0, nodes);
}

std::deque<LUT_Tree_Node*> LUT_Tree::get_level_nodes(int level){
    std::deque<LUT_Tree_Node*> nodes;
    get_level_nodes(level, nodes);
    return nodes;
}


// IO Operators
std::ostream& lut_ldpc::operator<<(std::ostream &os, const LUT_Tree &t){
    // Write tree header
    os << t.type << " " << t.num_leaves << std::endl;
    // Write Tree recursively
    t.root->serialize_tree(os);
    return os;
}

std::ostream& lut_ldpc::operator<<(std::ostream &os, const Array<Array<LUT_Tree>> &a){
    os << a.size() << std::endl;
    for(int ii=0; ii< a.size(); ii++){
        os << a(ii).size() << std::endl;
        for(int jj=0; jj<a(ii).size(); jj++){
            os << a(ii)(jj);
        }
    }
    return os;
}



std::istream& lut_ldpc::operator>>(std::istream &is, LUT_Tree &tree){
    
    std::string line;
    std::stringstream s;
    int t, numl;
    
    // Get Tree Type and number of leaves
    getline(is, line);
    s << line;
    s >> t >> numl;
    it_assert(!s.fail(),
              "std::istream& itpp::operator>>(std::istream &, LUT_Tree&): Error reading node type and  number of leaves");
    it_assert((t >= 0 && t < LUT_Tree::num_tree_types) && (numl >= 0),
              "std::istream& itpp::operator>>(std::istream &, LUT_Tree&): Wrong tree data");
    s.seekg(0, std::ios::end);
    s.clear();
    
    // Write to tree
    tree = LUT_Tree(); //Destroy a potentially existing tree
    tree.num_leaves = numl;
    tree.type = static_cast<LUT_Tree::tree_type_t>(t);
    tree.root = LUT_Tree_Node::deserialize_tree(is);
    return is;
}

std::istream& lut_ldpc::operator>>(std::istream &is, Array<Array<LUT_Tree>> &a){
    std::string line;
    std::stringstream s;
    int deg, num_trees_iter;
    
    // Get number of trees
    getline(is, line);
    s << line;
    s >> num_trees_iter;
    it_assert(!s.fail(),
              "std::istream& itpp::operator>>(std::istream &, Array<Array<LUT_Tree>>&): Error reading number of trees in dimension 1");
    it_assert((num_trees_iter >= 0) ,
              "std::istream& itpp::operator>>(std::istream &, Array<Array<LUT_Tree>>&): Wrong tree data");
    s.seekg(0, std::ios::end);
    s.clear();
    
    // Initialize a new array
    a = Array<Array<LUT_Tree>>(num_trees_iter);
    for(int ii=0; ii< num_trees_iter; ii++){
        getline(is, line);
        s << line;
        s >> deg;
        it_assert(!s.fail(),
                  "std::istream& itpp::operator>>(std::istream &, Array<Array<LUT_Tree>>&): Error reading number of trees in dimension 2");
        it_assert((deg >= 1) ,
                  "std::istream& itpp::operator>>(std::istream &, Array<Array<LUT_Tree>>&): Wrong tree data");
        s.seekg(0, std::ios::end);
        s.clear();
        a(ii) = Array<LUT_Tree>(deg);
        for(int jj=0; jj<deg; jj++){
            is >> a(ii)(jj) ;
        }
    }
    return is;
}

