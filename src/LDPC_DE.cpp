/*!
 * \file LDPC_DE.cpp
 * \brief Density Evolution for discrete LUT message passing decoding of LDPC codes
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

#include "LDPC_DE.hpp"
#include "TextTable.hpp"


using namespace lut_ldpc;
using namespace itpp;
using namespace std;


// ----------------------------------------------------------------------
// LDPC_DE
// ----------------------------------------------------------------------

void LDPC_DE::set_bisec_window(double tmin, double tmax){
    thr_min = tmin;
    thr_max = tmax;
    return;
}


int LDPC_DE::bisec_search(double & thr){
    it_assert(init_flag, "LDPC_DE::bisec_search(): Object not initialized");
    int ach = -1;
    bool converged = false;
    int ii=0;
    double sig = -1.0;
    double thr_min_ = thr_min;
    double thr_max_ = thr_max;
    while (converged==false && ii<maxiter_bisec) {
        // Update sig and check if achievable
        switch(mean_mode) {
            case ARI:
                sig = (thr_max_ + thr_min_)/2;
                break;
            case GEO:
                sig = std::sqrt(thr_max_*thr_min_);
                break;
            default:
                it_error("LDPC_DE::bisec_search(): Mean mode not defined");
                
        }
        it_info_debug("Evolving with threshold " << sig << "..." << std::endl);
        ach = evolve(sig);
        it_info_debug("Done, exit code = " << ach << std::endl);
        
        // If achievable and within precision, declare convergence
        if( (thr_max_-thr_min_ < thr_prec) && ach>=0 )
            converged = true;
        // If the bisection window size is within the precision and we are still not converged, exit
     //   if( (thr_max_-thr_min_ < thr_prec) && ach<0 )
     //       return -1;
        // Set new boundaries depending on ach
        if(ach>=0)
            thr_min_ = sig;
        else
            thr_max_ = sig;
        
        ii++;
    }
    if(converged){
        thr = sig;
        return ii;
    }
    else{
        thr = 0;
        return -1;
    }
}




// ----------------------------------------------------------------------
// LDPC_DE_LUT
// ----------------------------------------------------------------------

LDPC_DE_LUT::LDPC_DE_LUT(const LDPC_Ensemble ens_,
                               int Nq_Cha_,
                               ivec Nq_Msg_vec_,
                               int maxiter_de_,
                               Array<Array<LUT_Tree>> var_tree_templates_,
                               Array<Array<LUT_Tree>> chk_tree_templates_,
                               const bvec& reuse_vec_ ,
                               double thr_prec_ ,
                               double Pe_max_,
                               int mean_mode_,
                               int maxiter_bisec_,
                               double LLR_max_ ,
                               int Nq_fine_,
                               const std::string& irregular_design_strategy_)
{
    ens = ens_;
    Nq_Cha = Nq_Cha_;
    Nq_Msg_vec = Nq_Msg_vec_;
    maxiter_de = maxiter_de_;
    
    if(length(reuse_vec_)>0)
        reuse_vec = reuse_vec_;
    else
        reuse_vec = zeros_b(maxiter_de);
    
    thr_prec = thr_prec_;
    Pe_max = Pe_max_;
    mean_mode = mean_mode_;
    maxiter_bisec = maxiter_bisec_;
    LLR_max = LLR_max_;
    Nq_fine = Nq_fine_;
    var_tree_templates = var_tree_templates_;
    
    if(chk_tree_templates_.size()>0){
        chk_tree_templates = chk_tree_templates_;
        min_lut = false;
    }
    else{
        min_lut = true;
    }
    
    max_ni_de_iters = 1;
    
    thr_max = rate_to_shannon_thr( ens.get_rate() );
    thr_min = thr_max* 1e-4;
    
    if( irregular_design_strategy_ == "individual"){
        irregular_design_strategy = INDIVIDUAL;
    }
    else if (irregular_design_strategy_ == "joint_level"){
        irregular_design_strategy = JOINT_LEVEL;
    }
    else if (irregular_design_strategy_ == "joint_root"){
        irregular_design_strategy = JOINT_ROOT;
    }
    else{
        it_error("Irregular Design Strategy " << irregular_design_strategy_ << " unknown!");
    }
    // Input checking
    
    init_flag = true;
}

void LDPC_DE_LUT::set_ensemble(const LDPC_Ensemble& ens){
    // If the DE object is initialized, make sure the new ensemble matches the active degrees of the old ensemble
    if(init_flag){
        it_assert(
            ( sum(to_ivec(this->ens.get_dc_act() != ens.get_dc_act())) == 0 ) &&
            ( sum(to_ivec(this->ens.get_dv_act() != ens.get_dv_act())) == 0 ) ,
            "LDPC_DE_LUT::set_ensemble(): Active degree mismatch!"
        );
        this->ens = ens;
    }
    else{
        it_error("LDPC_DE_LUT::set_ensemble(): Object not initialized");
    }
}

int LDPC_DE_LUT::evolve(double thr, bool var, bool chk, mat& P, vec& p){
    Array<Array<LUT_Tree>> var_trees;
    Array<Array<LUT_Tree>> chk_trees;
    return evolve(thr, var, chk, P, p, false, var_trees, chk_trees);
}

int LDPC_DE_LUT::evolve(double thr){
    mat P;
    vec p;
    Array<Array<LUT_Tree>> var_trees;
    Array<Array<LUT_Tree>> chk_trees;
    return evolve(thr, false, false, P, p, false, var_trees, chk_trees);
}


int LDPC_DE_LUT::evolve(double thr, bool var_trace, bool chk_trace, mat& P, vec& p, bool save_luts, Array<Array<LUT_Tree>>& var_trees, Array<Array<LUT_Tree>>& chk_trees){
    it_assert_debug(reuse_vec(0)==bin(0), "LDPC_DE_LUT::evolve(): Reuse not possible for initial iteration");
    it_assert(!var_trace || !chk_trace, "LDPC_DE_LUT::evolve(): Choose either variable or check node tracing");
    
    // For the output after the last variable node update is binary. The original Nq_Msg_vec is restored at the end of this function.
    Nq_Msg_vec = concat(Nq_Msg_vec, 2);
    
    double Pe;
    double Pe_old = 1.0;
    int ni_iters = 0;
    
    // Number of trees to be stored when save_luts == true
    int num_qtrees = maxiter_de-sum(to_ivec(reuse_vec));

    // Set var2chk message and channel pmf
    set_channel_pmf(thr);
    
    // Get ensemble parameters
    int dv_act, dc_act;
    vec lam, rho;
    ivec degree_lam, degree_rho;
    dv_act = ens.sget_lam(lam, degree_lam);
    dc_act = ens.sget_rho(rho, degree_rho);
    
    // Stores trees that are going to be reused
    Array<LUT_Tree> var_trees_iter(dv_act);
    Array<LUT_Tree> chk_trees_iter(dc_act);
    
    // initialize probability tracing output if needed
    if(var_trace) {
        P.set_size(0, ens.get_dv_act());
        p.set_size(0);
    }
    
    if(chk_trace){
        P.set_size(0, ens.get_dc_act());
        p.set_size(0);
    }
    if(save_luts){
        var_trees.set_size(num_qtrees);
        if(min_lut){
            chk_trees.set_size(0);
            chk_trees_iter.set_size(0);
        }
        else{
            chk_trees.set_size(num_qtrees);
        }
    }
    
    
    int var_tree_idx = 0;
    int chk_tree_idx = 0;
    
    int max_iter;
    // We want to keep the actual LUTs including the decision tree
    if(save_luts)
        max_iter = maxiter_de;
    else
        max_iter = maxiter_de-1;
    // evolve
    for(int ii=0; ii<max_iter; ii++){
        // Convergence check
        Pe = sum(pmf_var2chk.left(Nq_Msg_vec(ii)/2) );
        it_info_debug("Iteration " << ii << " Pe= " << Pe);
        if(Pe < Pe_max && !save_luts)
            return ii;
        if(Pe <= Pe_old) Pe_old = Pe;
        else ni_iters++;
        
        if(ni_iters >= max_ni_de_iters && !save_luts)
            return -1;
        
        // ====== CHK update
        
        if(chk_trace){
            double p_elem;
            vec P_row(ens.get_dc_act());
            chk_update_irr(P_row,p_elem, ii, chk_trees_iter);
            P.append_row(P_row);
            p = concat(p, p_elem);
        }
        else{
            chk_update_irr(ii, chk_trees_iter);
        }
        
        // ====== VAR update
        if(var_trace){
            double p_elem;
            vec P_row(ens.get_dv_act());
            var_update_irr(P_row,p_elem, ii, var_trees_iter);
            P.append_row(P_row);
            p = concat(p, p_elem);
        }
        else{
            var_update_irr(ii, var_trees_iter);
        }
        
        if(save_luts && reuse_vec(ii)==false){
            var_trees(var_tree_idx) = var_trees_iter;
            var_tree_idx++;
            if(min_lut == false){
                chk_trees(chk_tree_idx) = chk_trees_iter;
                chk_tree_idx++;
            }
            
        }
        
    }
    
    //remove pmfs
    if(save_luts){
        for(int ii=0; ii<var_trees.size(); ii++){
            for (int dd=0; dd<var_trees(ii).size(); dd++){
                var_trees(ii)(dd).reset_pmfs();
            }
        }
        for(int ii=0; ii<chk_trees.size(); ii++){
            for (int dd=0; dd<chk_trees(ii).size(); dd++){
                chk_trees(ii)(dd).reset_pmfs();
            }
        }
    }
    
    Nq_Msg_vec.del(Nq_Msg_vec.length()-1);
    if(!save_luts)
        return -1;
    else
        return max_iter;
}

bvec LDPC_DE_LUT::evolve_adaptive_reuse(double thr, vec& Pe_trace, double rel_increase_max, double rel_decrease_min, int reuse_max){
     
    double Pe;
    double Pe_old = 1.0;
    int ni_iters = 0;
   
    // Save initial reuse vec and reset upon exit of fuction
    bvec reuse_vec_old = reuse_vec; 
    // Set var2chk message and channel pmf
    set_channel_pmf(thr);
    
    // Get ensemble parameters
    int dv_act = ens.get_dv_act();
    int dc_act = ens.get_dc_act();

    // Stores trees that are going to be reused
    Array<LUT_Tree> var_trees_iter(dv_act);
    Array<LUT_Tree> chk_trees_iter(dc_act);
    
    // initialize probability tracing output if needed
    Pe_trace.set_size(0);

    // Initialize Reuse vector
    reuse_vec = zeros_b(reuse_vec.length());   

    // evolve
    int ii = 0;
    int num_reuse = 0;
    for(ii=0; ii<maxiter_de-1; ii++){
        // Convergence check
        Pe = sum(pmf_var2chk.left(Nq_Msg_vec(ii)/2) );
        double rel_decrease = (Pe_old-Pe)/Pe_old;
        it_info_debug("Iteration " << ii << " Pe= " << Pe << ", rel_decrease = " << rel_decrease);
        if(Pe < Pe_max)
            break; 
        if(Pe <= Pe_old) Pe_old = Pe;
        else ni_iters++;
        
        if(ni_iters >= max_ni_de_iters)
            break; 
        if(ii!=0){
            reuse_vec(ii)=1;
        }
        // Save olde density and try an iteration with reuse
        vec pmf_var2chk_old = pmf_var2chk;
        chk_update_irr(ii, chk_trees_iter);
        var_update_irr(ii, var_trees_iter);
        double Pe_new_reuse = sum(pmf_var2chk.left(Nq_Msg_vec(ii)/2));
        double Pe_old = sum(pmf_var2chk_old.left(Nq_Msg_vec(ii)/2));
        double rel_increase = (Pe_new_reuse-Pe_old)/Pe_old;
        if( rel_increase > rel_increase_max || -rel_increase < rel_decrease_min || num_reuse > reuse_max){
            it_info_debug(" Relative decrease " << -rel_increase << " from Pe = " << Pe_old << " to Pe = " << Pe_new_reuse 
                    << ", reuses = " << num_reuse);
            reuse_vec(ii)=0;
            pmf_var2chk = pmf_var2chk_old;
            chk_update_irr(ii, chk_trees_iter);
            var_update_irr(ii, var_trees_iter);
            num_reuse = 0;
        }
        else{
            num_reuse++;
        }
    }
    
    bvec reuse_vec_return = reuse_vec.left(ii);
    reuse_vec = reuse_vec_old;
    return reuse_vec_return;
}

double LDPC_DE_LUT::get_lam2stable(double sig){
    return get_lam2stable_lut(sig, this->ens.get_chk_degree_dist(), Nq_Cha,  Nq_Msg_vec(0));  
}

void LDPC_DE_LUT::set_channel_pmf(double sig){
    double delta = 2*LLR_max/Nq_fine;
    //double delta = 2*(2/sqr(sig) + 4*2/sig)/Nq_fine;
    vec pmf_channel_fine = get_gaussian_pmf(2/sqr(sig), 2/sig, Nq_fine, delta);
    
    ivec Q_out;
    
    (void) quant_mi_sym(pmf_cha,     Q_out, pmf_channel_fine, Nq_Cha, true);
    (void) quant_mi_sym(pmf_var2chk, Q_out, pmf_channel_fine, Nq_Msg_vec(0), true);
    
    
    return;
}

void LDPC_DE_LUT::chk_update_irr(vec& P_row, double& Pe, int iter,  Array<LUT_Tree>& prev_trees_chk){
    
    int dc_act = ens.get_dc_act();
    vec rho = ens.sget_rho();
    ivec degree_rho = ens.sget_degree_rho();
    
    pmf_chk2var = zeros(Nq_Msg_vec(iter));
    
    P_row.set_size(dc_act);
    Pe = 0;
    
    if(min_lut){
        for(int dd=0; dd<dc_act; dd++){
            vec p_tmp = chk_update_minsum(pmf_var2chk, degree_rho(dd));
            P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
            Pe += rho(dd)*P_row(dd);
            pmf_chk2var =  pmf_chk2var + rho(dd)*p_tmp;
        }
    }
    else{
        if(reuse_vec(iter)){
            for(int dd=0; dd<dc_act; dd++){
                prev_trees_chk(dd).set_leaves(pmf_var2chk, pmf_cha);
                vec p_tmp = prev_trees_chk(dd).update(true);
                P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
                Pe += rho(dd)*P_row(dd);
                pmf_chk2var = pmf_chk2var + rho(dd)*p_tmp;
            }
        }
        else{
            // Copy trees from templates and prepare for updates
            for(int dd=0; dd<dc_act; dd++){
                LUT_Tree tree = chk_tree_templates(iter)(dd);
                tree.set_leaves(pmf_var2chk, pmf_cha);
                tree.set_resolution(Nq_Msg_vec(iter), Nq_Msg_vec(iter), Nq_Cha);
                prev_trees_chk(dd) = tree;
            }
            
            // Update the tree LUTs
            switch (irregular_design_strategy) {
                case INDIVIDUAL:
                    for(int dd=0; dd<dc_act; dd++){
                        vec p_tmp = prev_trees_chk(dd).update();
                        P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
                        Pe += rho(dd)*P_row(dd);
                        pmf_chk2var = pmf_chk2var + rho(dd)*p_tmp;
                    }
                    break;
                case JOINT_LEVEL:
                    pmf_chk2var = joint_level_irr_lut_design(rho, degree_rho, prev_trees_chk, P_row, Pe);
                    pmf_chk2var.zeros();
                    for(int dd=0; dd<dc_act; dd++){
                        vec p_tmp = prev_trees_chk(dd).update(true);
                        P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
                        Pe += rho(dd)*P_row(dd);
                        pmf_chk2var = pmf_chk2var + rho(dd)*p_tmp;
                    }
                    break;
                case JOINT_ROOT:
                    pmf_chk2var = joint_root_irr_lut_design(rho, degree_rho, prev_trees_chk, P_row, Pe);
                    pmf_chk2var.zeros();
                    for(int dd=0; dd<dc_act; dd++){
                        vec p_tmp = prev_trees_chk(dd).update(true);
                        P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
                        Pe += rho(dd)*P_row(dd);
                        pmf_chk2var = pmf_chk2var + rho(dd)*p_tmp;
                    }
                    break;
                default:
                    it_error("LDPC_DE_LUT::chk_update_irr(): Irregular design strategy not known!");
                    break;
            }
            
        }
    }
}




void LDPC_DE_LUT::var_update_irr(vec& P_row, double& Pe, int iter,  Array<LUT_Tree>& prev_trees_var){
    
    int dv_act = ens.get_dv_act();
    vec lam = ens.sget_lam();
    ivec degree_lam = ens.sget_degree_lam();
    
    pmf_var2chk = zeros(Nq_Msg_vec(iter+1));
    
    P_row.set_size(dv_act);
    Pe = 0;
    
    if(reuse_vec(iter)){
        for(int dd=0; dd<dv_act; dd++){
            prev_trees_var(dd).set_leaves(pmf_chk2var, pmf_cha);
            vec p_tmp = prev_trees_var(dd).update(true);
            P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
            Pe += lam(dd)*P_row(dd);
            pmf_var2chk = pmf_var2chk + lam(dd)* p_tmp;
        }
    }
    else{
        // Copy trees from templates and prepare for updates
        for(int dd=0; dd<dv_act; dd++){
            LUT_Tree tree = var_tree_templates(iter)(dd);
            tree.set_leaves(pmf_chk2var, pmf_cha);
            tree.set_resolution(Nq_Msg_vec(iter), Nq_Msg_vec(iter+1), Nq_Cha);
            prev_trees_var(dd) = tree;
        }
        // Update the tree LUTs
        switch (irregular_design_strategy) {
            case INDIVIDUAL:
                for(int dd=0; dd<dv_act; dd++){
                    vec p_tmp = prev_trees_var(dd).update();
                    P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
                    Pe += lam(dd)*P_row(dd);
                    pmf_var2chk = pmf_var2chk + lam(dd)* p_tmp;
                }
                break;
            case JOINT_LEVEL:
                pmf_var2chk = joint_level_irr_lut_design(lam, degree_lam, prev_trees_var, P_row, Pe);
                pmf_var2chk.zeros();
                for(int dd=0; dd<dv_act; dd++){
                    vec p_tmp = prev_trees_var(dd).update(true);
                    P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
                    Pe += lam(dd)*P_row(dd);
                    pmf_var2chk = pmf_var2chk + lam(dd)* p_tmp;
                }
                break;
            case JOINT_ROOT:
                pmf_var2chk = joint_root_irr_lut_design(lam, degree_lam, prev_trees_var, P_row, Pe);
                pmf_var2chk.zeros();
                for(int dd=0; dd<dv_act; dd++){
                    vec p_tmp = prev_trees_var(dd).update(true);
                    P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
                    Pe += lam(dd)*P_row(dd);
                    pmf_var2chk = pmf_var2chk + lam(dd)* p_tmp;
                }
                break;
            default:
                it_error("LDPC_DE_LUT::var_update_irr(): Irregular design strategy not known!");
                break;
        }
        
    }
}


void LDPC_DE_LUT::get_quant_bound(double sig, vec& qb_Cha, vec& qb_Msg) const{
    double delta = 2*LLR_max/Nq_fine;
    //double delta = 2*(2/sqr(sig) + 4*2/sig)/Nq_fine;
    vec pmf_channel_fine = get_gaussian_pmf(2/sqr(sig), 2/sig, Nq_fine, delta);
    int M = Nq_fine;
    int K;
    
    ivec Q_out;
    vec p;
    
    // calculate quantization boundaries for channel LLRs
    (void) lut_ldpc::quant_mi_sym(p,     Q_out, pmf_channel_fine, Nq_Cha, true);
    K = Nq_Cha;
    Q_out = (Q_out.right(M/2)-K/2);
    qb_Cha.set_size(K/2-1);
    int label = 0;
    for(int mm=0; mm<M/2; mm++){
        if(Q_out(mm) > label ){
            qb_Cha(label) = mm*delta;
            label++;
        }
    }
    qb_Cha = concat(concat( -fliplr(qb_Cha), to_vec(0)), qb_Cha);
    
    // calculate quantization boundaries for channel LLRs
    (void) lut_ldpc::quant_mi_sym(p,     Q_out, pmf_channel_fine, Nq_Msg_vec(0), true);
    K = Nq_Msg_vec(0);
    Q_out = (Q_out.right(M/2)-K/2);
    qb_Msg.set_size(K/2-1);
    label = 0;
    for(int mm=0; mm<M/2; mm++){
        if(Q_out(mm) > label ){
            qb_Msg(label) = mm*delta;
            label++;
        }
    }
    qb_Msg = concat(concat( -fliplr(qb_Msg), to_vec(0)), qb_Msg);

    return;
    
}





void LDPC_DE_LUT::get_lut_trees(Array<Array<LUT_Tree>>& var_trees, Array<Array<LUT_Tree>>& chk_trees, double sig){
    it_assert_debug(reuse_vec(0)==bin(0), "LDPC_DE_MinLUT::get_var_lut_trees(): Reuse not possible for initial iteration");
    mat P;
    vec p;
    evolve(sig, false, false, P, p, true, var_trees, chk_trees);
}

// ----------------------------------------------------------------------
// LDPC_DE_BP
// ----------------------------------------------------------------------
LDPC_DE_BP::LDPC_DE_BP(LDPC_Ensemble ens_,   int Nb_, double Lmax_){
    Nb = Nb_;
    Lmax = Lmax_;
    N = pow2i(Nb-1);
    delta = 2*Lmax/(2*N+1);
    ens = ens_;
    mean_mode = LDPC_DE::ARI;
    Nfft = pow2i(1+std::ceil(std::log2(2*N+1)));
    
    // Space allocation
    pmf_LLR.set_size(2*N+2);
    pmf_chk2var.set_size(2*N+2);
    pmf_var2chk.set_size(2*N+2);
    
    idx_support_sym = linspace_fixed_step(-N, N);
    idx_support = idx_support_sym + N;
    
    support = concat(idx_support_sym*delta, std::numeric_limits<double>::max());
    var_conv_weight = exp(-.5* idx_support_sym*delta);
    
    
    
    // Exit conditions
    max_ni_de_iters = 1;
    maxiter_de = 1000;
    Pe_max = 1e-9;
    maxiter_bisec = 50;
    thr_prec = 1e-4;
    
    
    // Ensemble specific configurations
    ens = ens_;
    thr_max = rate_to_shannon_thr(ens.get_rate());
    thr_min = thr_max/1e3;
    
    
    set_tq_tables();
    
    
    setup_flag = true;
    init_flag = true;
}

void LDPC_DE_BP::set_ensemble(const LDPC_Ensemble& ens){
        this->ens = ens;
}


int LDPC_DE_BP::evolve(double thr){
    mat P;
    vec p;
    return  evolve(thr, false, false, P, p);
}

int LDPC_DE_BP::evolve(double thr,  bool var_trace, bool chk_trace, mat& P, vec& p){
    it_assert(init_flag & setup_flag, "LDPC_DE_BP::evolve(): Object has not been initialized");
    it_assert(!var_trace || !chk_trace, "LDPC_DE_BP::evolve(): Choose either variable or check node tracing");
    int ni_iters = 0; // number of DE iterations where Pe did not improve
    double Pe_old = 1.0;
    double Pe = 1.0;
    
    // reset densities
    pmf_LLR = concat(get_gaussian_pmf(2/sqr(thr), 2/thr, 2*N+1, delta), 0.0);
    pmf_var2chk    = pmf_LLR;
    
    // initialize probability tracing output if needed
    if(var_trace) {
            P.set_size(0, ens.get_dv_act());
            p.set_size(0);
    }

    if(chk_trace){
        P.set_size(0, ens.get_dc_act());
        p.set_size(0);
    }
    
    for(int ii=1; ii<maxiter_de; ii++){
        if(chk_trace){
            double p_elem;
            vec P_row(ens.get_dc_act());
            chk_update_irr(true, P_row,p_elem);
            P.append_row(P_row);
            p = concat(p, p_elem);
        }
        else{
            chk_update_irr();
        }

        
        //=== Step 2: VN update
        if (var_trace) {
                double p_elem;
                vec P_row(ens.get_dv_act());
                var_update_irr(var_trace, P_row,p_elem);
                P.append_row(P_row);
                p = concat(p, p_elem);
        }
        else{
            var_update_irr();
        }

        
        
        Pe = sum(pmf_var2chk.get(0, N-1))+ .5*pmf_var2chk(N);
        //=== convergence check
        it_info_debug("Iteration "<< ii <<  " Pe = " << Pe);
        if( Pe  < Pe_max )
            return ii+1;
        
        if(Pe < Pe_old) Pe_old = Pe;
        else ni_iters++;
        
        if(ni_iters >= max_ni_de_iters)
            return -1;
        
    }
    return -1;
}

double LDPC_DE_BP::get_lam2stable(double sig){
    return get_lam2stable_cbp(sig, this->ens.get_chk_degree_dist());  
}

vec LDPC_DE_BP::pmf_plus(const vec& pmf){
    it_assert(length(pmf)== 2*N+2, "SMLDPC_DE_PMF::pmf_plus(): Input Dimension does not match object resolution");
    // the positive part includes the mass at zero and infinity
    vec v(N+2);
    v(0) = pmf(N);
    for(int nn=1; nn<N+1; nn++){
        v(nn) = pmf(N+nn) + pmf(N-nn);
    }
    v(N+1) = pmf(2*N+1);
    
    return v;
}

vec LDPC_DE_BP::pmf_minus(const vec& pmf){
    it_assert(length(pmf)== 2*N+2, "SMLDPC_DE_PMF::pmf_plus(): Input Dimension does not match object resolution");
    vec v(N+2);
    v(0)=0;
    for(int nn=1; nn<N+1; nn++){
        v(nn) = pmf(N+nn) - pmf(N-nn);
    }
    v(N+1) = pmf(2*N+1);
    return v;
}

vec LDPC_DE_BP::pmf_orig(const vec& pmf_p, const vec& pmf_m){
    it_assert(length(pmf_p)== N+2 && length(pmf_m)== N+2, "SMLDPC_DE_PMF::pmf_plus(): Input Dimension does not match object resolution");
    vec v(2*N+2);
    // reconstruct negative part
    for(int nn=1; nn<N+1; nn++){
        v(N-nn) = .5* (pmf_p(nn) - pmf_m(nn));
    }
    // mass at zero
    v(N) = pmf_p(0);
    // reconstruct positive part (including infinity mass)
    for(int nn=1; nn<N+2; nn++){
        v(N+nn) = .5* (pmf_p(nn) + pmf_m(nn));
    }
    //v(2*N+1) = pmf_p(N+1);
    return v;
}



void LDPC_DE_BP::chk_update_irr(bool trace, vec& P_row, double& Pe){
    vec pmf_out_tmp;
    vec pmf_out_tmp_p;
    vec pmf_out_tmp_m;
    
    
    pmf_chk2var.zeros();
    
    vec pmf_in_p = pmf_plus(pmf_var2chk);
    vec pmf_in_m = pmf_minus(pmf_var2chk);
    
    // Get and sort check node distribution
    ivec dc_vec = ens.sget_degree_rho();
    int dc_act = ens.get_dc_act();
    vec chk_degree_dist = ens.sget_rho();
    
    ivec idx = sort_index(dc_vec);
    dc_vec = dc_vec.get(idx);
    chk_degree_dist = chk_degree_dist.get(idx);
    
    // Set initial output pmf to input pmf. The output pmfs will be concolved against the input pmf in the following
    pmf_out_tmp_p = pmf_in_p;
    pmf_out_tmp_m = pmf_in_m;
    
    int dc_tmp = 2; // determines the check node degree currently represented by the *_tmp distribution
    
    for(int jj=0; jj<dc_act; jj++){
        if(chk_degree_dist(jj)!= 0 ||  trace){
            if(dc_vec(jj)==1 || dc_vec(jj)==2){
                // Output is equal to input as set before => do nothing
            }
            else{
                int conv_times =  dc_vec(jj)-dc_tmp;
                it_assert(conv_times >0, "SMLDPC_DE_PMF::chk_update_irr(): Check node degrees must be strictly increasing");
                for(int ii=0; ii< conv_times; ii++){
                    chk_update_convolve(pmf_in_p, pmf_in_m, pmf_out_tmp_p, pmf_out_tmp_m);
                    dc_tmp++;
                }
            }
            pmf_out_tmp = pmf_orig(pmf_out_tmp_p, pmf_out_tmp_m);
            pmf_chk2var += chk_degree_dist(jj)*pmf_out_tmp;
            if(trace) P_row(jj) = sum(pmf_out_tmp.left(N)) + .5*pmf_out_tmp.get(N);
        }
    }
    if(trace) Pe = sum(pmf_chk2var.left(N)) + .5*pmf_chk2var.get(N);
    return;
}


void LDPC_DE_BP::var_update_irr(bool trace, vec& P_row, double& Pe){
    
    pmf_var2chk.zeros();
    
    
    //====  Get and sort variable node distribution
    ivec dv_vec = ens.sget_degree_lam();
    vec var_degree_dist = ens.sget_lam();
    int dv_act = ens.get_dv_act();
    
    ivec idx = sort_index(dv_vec);
    dv_vec = dv_vec.get(idx);
    var_degree_dist = var_degree_dist.get(idx);
    
    vec pmf_tmp = pmf_LLR;
    int dv_tmp = 1; // determines the variable node degree currently represented by the pmf_tmp distribution
    
    
    
    for(int jj=0; jj<dv_act; jj++){
        if(var_degree_dist(jj)!= 0 ||  trace){
            if(dv_vec(jj)==1){
                // Output is equal to input as set before => do nothing
            }
            else{
                int conv_times =  dv_vec(jj)-dv_tmp;
                it_assert(conv_times >0, "LDPC_DE_BP::chk_update_irr(): Variable node degrees must be strictly increasing");
                for(int ii=0; ii< conv_times; ii++){
                    var_update_convolve(pmf_chk2var, pmf_tmp);
                    dv_tmp++;
                }
            }


            // upate distribution
            pmf_var2chk = pmf_var2chk + var_degree_dist(jj)*pmf_tmp;

            if(trace) {
                P_row(jj) = sum(pmf_tmp.get(0, N-1)) + .5*pmf_tmp.get(N);
            }
        }
    }
    if(trace) Pe = sum(pmf_var2chk.left(N)) + .5*pmf_var2chk.get(N);
    return;
}


vec LDPC_DE_BP::fft_preprocess(const vec& x) const {
    vec y = zeros(Nfft);
    y.set_subvector(0, x.get(N, 2*N));
    y.set_subvector(Nfft-N, x.get(0,N-1));
    return y;
}

vec LDPC_DE_BP::ifft_postprocess(const vec& x) const {
    vec y = concat(x.right(N), x.left(N+1));
    y = concat(y, 1-sum(y));
    return y;
}

void LDPC_DE_BP::numeric_postprocessing(vec& x) const {
    for(int ii=0; ii<x.length(); ii++) if(x(ii)<0)  x(ii)=0;
    x = x/sum(x);
}



void LDPC_DE_BP::var_update_convolve(const vec &pmf_in, vec &pmf_out){
//    //############# general version without weighting
//    it_assert(Nfft > 2*N+2, "SMLDPC_DE_PMF::var_update_convolve(): FFT size too small");
//    // Scale and make even pmfs
//    vec a_hat =  pmf_in.left(2*N+1);
//    vec b_hat =  pmf_out.left(2*N+1);
//    // convolve in frequency domain
//    vec pmf_tmp = ifft_real( elem_mult( fft_real(a_hat, Nfft), fft_real(b_hat,Nfft) ) );
//    // Cut out part correspondin to input pmf's support, neglecting the probability mass in [-2*N,-(N+1)] We obtain mass in [N+1, 2N] later
//    pmf_out = pmf_tmp.get(N, 3*N);
//    pmf_out(0) += sum(pmf_tmp.get(0, N-1));
//    pmf_out(2*N) += sum(pmf_tmp.get(3*N+1, 4*N+1));
//    // Compute mass at inf
//    pmf_out = concat(pmf_out, 1-sum(pmf_out));
    
    //############# symmetric optimistic version with weighting
    it_assert(Nfft > 2*N+2, "SMLDPC_DE_PMF::var_update_convolve(): FFT size too small");
    // Scale and make even pmfs
    vec a_hat =  elem_mult(var_conv_weight, pmf_in.left(2*N+1));
    vec b_hat =  elem_mult(var_conv_weight, pmf_out.left(2*N+1));
    // convolve in frequency domain
    vec pmf_tmp = ifft_real( elem_mult( fft_real(a_hat, Nfft), fft_real(b_hat,Nfft) ) );
   
    // Cut out part correspondin to input pmf's support, neglecting the probability mass in [-2*N,-(N+1)] We obtain mass in [N+1, 2N] later
    pmf_tmp = pmf_tmp.get(N, 3*N);
    // Undo scaling
    pmf_out = elem_div(pmf_tmp, var_conv_weight);
    // Compute mass at inf
    pmf_out = concat(pmf_out, 1- sum(pmf_out));
    
//    //  numeric_postprocessing(pmf_out);
    
//    
//    it_assert(Nfft > 4*N, "SMLDPC_DE_PMF::var_update_convolve(): FFT size too small");
//    // Scale and make even pmfs
//    double a_inf = pmf_in(2*N+1);
//    double b_inf = pmf_out(2*N+1);
//
//    // convolve in frequency domain
//    vec pmf_tmp = ifft_real( elem_mult( fft_real(pmf_in.left(2*N+1), Nfft), fft_real(pmf_out.left(2*N+1),Nfft) ) );
//    
//
//    
//    // Cut out part correspondin to input pmf's support, neglecting the probability mass in [-2*N,-(N+1)] We obtain mass in [N+1, 2N] later
//    pmf_out = pmf_tmp.get(N, 3*N);
//    // Add Masses at +-N
//    pmf_out(0) = pmf_out(0) + sum(pmf_tmp.get(0,N-1));
//    pmf_out(2*N) = pmf_out(2*N) + sum(pmf_tmp.get(3*N+1,Nfft-1));
//    // Compute mass at inf
//    pmf_out = concat(pmf_out, a_inf + b_inf - a_inf*b_inf);
//  //  std::cout << sum(pmf_out) << std::endl;
//
//    numeric_postprocessing(pmf_out);
//    

    
    return;
}

void LDPC_DE_BP::chk_update_convolve(const vec &pmf_in_p, const vec &pmf_in_m, vec &pmf_out_p, vec &pmf_out_m){
    
    
    vec a_plus_fin =  pmf_in_p.left(N+1);
    vec a_minus_fin = pmf_in_m.left(N+1);
    vec b_plus_fin =  pmf_out_p.left(N+1);
    vec b_minus_fin = pmf_out_m.left(N+1);
    
    vec Aplus =  concat( concat(to_vec(sum(a_plus_fin)),  sum(a_plus_fin)  - cumsum(a_plus_fin.left(N)))  + pmf_in_p.right(1).get(0),  to_vec(0));
    vec Aminus = concat( concat(to_vec(sum(a_minus_fin)), sum(a_minus_fin) - cumsum(a_minus_fin.left(N))) + pmf_in_m.right(1).get(0),  to_vec(0));
    vec Bplus =  concat( concat(to_vec(sum(b_plus_fin)),  sum(b_plus_fin)  - cumsum(b_plus_fin.left(N)))  + pmf_out_p.right(1).get(0), to_vec(0));
    vec Bminus = concat( concat(to_vec(sum(b_minus_fin)), sum(b_minus_fin) - cumsum(b_minus_fin.left(N))) + pmf_out_m.right(1).get(0), to_vec(0));
    
    vec c_plus   = zeros(N+2);
    vec c_minus  = zeros(N+2);
    for (int i=0; i<=N; i++) {
        for (int k=0; k<=K; k++) {
            
            if ( (i-k)>=0 )
            {
                c_plus(i-k) += pmf_in_p(i)*(Bplus(tq(i,k+1)) - Bplus(tq(i,k)))
                + pmf_out_p(i)*(Aplus(tq2(i,k+1)) - Aplus(tq2(i,k)));
                
                c_minus(i-k) += pmf_in_m(i)*(Bminus(tq(i,k+1)) - Bminus(tq(i,k)))
                + pmf_out_m(i)*(Aminus(tq2(i,k+1)) - Aminus(tq2(i,k)));
            }
        }
    }

    c_plus(N+1) = pmf_in_p(N+1)*pmf_out_p(N+1);
    c_minus(N+1) = pmf_in_p(N+1)*pmf_out_p(N+1);
    
    pmf_out_p = c_plus;
    pmf_out_m = c_minus;
    
    
    return;
}



imat LDPC_DE_BP::gen_Q_table(){
    imat Q(N+1,N+1);
    for(int ii=0; ii<N+1; ii++){
        for(int jj=0; jj<N+1; jj++){
            Q(ii,jj) = floor_i( 2*std::atanh( std::tanh(.5*ii*delta) * std::tanh(.5*jj*delta) )/delta +.5);
        }
    }
    // add a column and row for infinity inputs
    Q.append_col(linspace_fixed_step(0, N));
    Q.append_row(linspace_fixed_step(0, N+1));
    Q(N+1,N+1) = N+1;
    return Q;
}

void LDPC_DE_BP::set_tq_tables(){
    // calculate Q table
    imat Q = gen_Q_table();
    //update table size
    K = ceil_i(std::log(2)/delta-.5);
    tq.set_size(N+1, K+2);
    tq2.set_size(N+1, K+2);
    for(int ii=0; ii<N+1; ii++){
        for(int kk=0; kk<K+2; kk++){
            if(kk==0){
                tq(ii,kk) = N+1;
                tq2(ii,kk) = N+1;
            }
            else{
                ivec row_ii= Q.get_row(ii);
                ivec q_indices = find(row_ii >= (ii -(kk-1) ) );
                int q_idx;
                if(q_indices.length()>0)
                    q_idx = q_indices(0);
                else
                    q_idx = std::numeric_limits<int>::min();
                
                tq(ii,kk)  = std::max(ii,   q_idx);
                tq2(ii,kk) = std::max(ii+1, q_idx);
            }
        }
    }
    return;
    
}

// ----------------------------------------------------------------------
// Other functions
// ----------------------------------------------------------------------














vec lut_ldpc::chk_update_minsum(const vec& p_in, int dc){
    
    int N = length(p_in);
    it_assert_debug(mod(N, 2)==0, "chk_update_minsum(): Input pmf must have even length");
    vec c_plus(N/2);
    vec c_minus(N/2);
    
    
    vec a_plus  = pmf_plus(p_in);
    vec a_minus = pmf_minus(p_in);
    vec b_plus  = a_plus;
    vec b_minus = a_minus;
    
    for(int dd=1; dd<dc-1; dd++){
        c_plus.clear();
        c_minus.clear();
        for(int ii=0; ii<N/2; ii++){
            for(int jj=0; jj<N/2; jj++){
                int kk = std::min(ii,jj);
                c_plus(kk) += a_plus(ii)*b_plus(jj);
                c_minus(kk) += a_minus(ii)*b_minus(jj);
            }
        }
        a_plus = c_plus;
        a_minus = c_minus;
    }
    
    return pmf_join(c_plus, c_minus);
}

inline vec lut_ldpc::pmf_plus(const vec& pmf){
    int N = pmf.length();
    it_assert_debug(mod(N, 2)==0, "pmf_plus(): Input pmf must have even length");
    vec pmf_p(N/2);
    for(int nn=0; nn<N/2; nn++){
        pmf_p(nn) = pmf(N/2+nn) + pmf(N/2-1-nn);
    }
    return pmf_p;
}

inline vec lut_ldpc::pmf_minus(const vec& pmf){
    int N = pmf.length();
    it_assert_debug(mod(N, 2)==0, "pmf_minus(): Input pmf must have even length");
    vec pmf_m(N/2);
    for(int nn=0; nn<N/2; nn++){
        pmf_m(nn) = pmf(N/2+nn) - pmf(N/2-1-nn);
    }
    return pmf_m;
}

inline vec lut_ldpc::pmf_join(const vec& pmf_p, const vec& pmf_m){
    int N = length(pmf_p);
    it_assert_debug(N == length(pmf_m), "pmf_join(): Length of input does not match");
    N = N*2;
    vec p_out(N);
    for(int nn=0; nn<N/2; nn++){
        p_out(N/2+nn)   = .5*( pmf_p(nn) + pmf_m(nn) );
        p_out(N/2-1-nn) = .5*( pmf_p(nn) - pmf_m(nn) );
    }
    return p_out;
}


void lut_ldpc::get_lut_tree_templates(const std::string& tree_method, const LDPC_Ensemble& ens, ivec Nq_Msg, int Nq_Cha, bool minLUT, Array<Array<LUT_Tree> >& var_luts, Array<Array<LUT_Tree> >& chk_luts )
{
    
    int max_iters = length(Nq_Msg);
    std::stringstream ss(tree_method);
    std::string tm;
    std::string filename;
    
    std::getline(ss, tm, '=');
    std::getline(ss, filename);
    ss.str(std::string());
    ss.clear();
    
    
    
    
    ivec var_deg = ens.sget_degree_lam();
    int dv_act  = ens.get_dv_act();
    ivec chk_deg = ens.sget_degree_rho();
    int dc_act  = ens.get_dc_act();
    
    
    if(tm=="filename"){
        it_assert(filename.length()>0, "lut_ldpc::get_lut_tree_templates(): To load the treestructure from a file, specify the tree_method input as 'filename=<name of file>'");
        boost::property_tree::ptree param_tree;
        boost::property_tree::ini_parser::read_ini(filename, param_tree);
        boost::optional<boost::property_tree::ptree&> cur_tree; // stores the current tree structure
        
        // Load trees for the initial iteration
        cur_tree   = param_tree.get_child_optional( "var_iter_000" );
        it_assert(cur_tree, "lut_ldpc::get_lut_tree_templates(): Error reading variable node tree from file!");
        
        var_luts.set_length(max_iters);
        var_luts(0).set_length(dv_act);
        for(int dd=0; dd<dv_act; dd++){
            ss << "var_deg_" << std::setfill('0') << std::setw(3) << var_deg(dd);
            boost::optional<std::string> cur_string = cur_tree->get_optional<std::string>( ss.str() );
            it_assert(cur_string, "lut_ldpc::get_lut_tree_templates(): Could not load tree string for variable node degree " << var_deg(dd) << " at iteration 0");
            var_luts(0)(dd) = LUT_Tree(*cur_string, LUT_Tree::VARTREE);
            it_assert(var_luts(0)(dd).get_num_leaves() == var_deg(dd), "lut_ldpc::get_lut_tree_templates(): LUT tree does not match node degree!");
            var_luts(0)(dd).set_resolution(Nq_Msg(0), Nq_Msg(1), Nq_Cha);
            ss.str(std::string());
            ss.clear();
        }
        for(int ii=1; ii< max_iters-1; ii++){
            var_luts(ii).set_length(dv_act);
            ss << "var_iter_" << std::setfill('0') << std::setw(3) << ii;
            cur_tree = param_tree.get_child_optional( ss.str() );
            ss.str(std::string());
            ss.clear();
            if(cur_tree){
                for(int dd=0; dd<dv_act; dd++){
                    ss << "var_deg_" << std::setfill('0') << std::setw(3) << var_deg(dd);
                    boost::optional<std::string> cur_string  = cur_tree->get_optional<std::string>( ss.str() );
                    it_assert(cur_string, "lut_ldpc::get_lut_tree_templates(): Could not load tree string for variable node degree " << var_deg(dd) << " at iteration " << ii << ".");
                    var_luts(ii)(dd) = LUT_Tree(*cur_string, LUT_Tree::VARTREE);
                    it_assert(var_luts(ii)(dd).get_num_leaves() == var_deg(dd), "lut_ldpc::get_lut_tree_templates(): LUT tree does not match node degree!");
                    var_luts(ii)(dd).set_resolution(Nq_Msg(ii), Nq_Msg(ii+1), Nq_Cha);
                    ss.str(std::string());
                    ss.clear();
                }
            }
            else{
                for(int dd=0; dd<dv_act; dd++){
                    var_luts(ii)(dd) = var_luts(ii-1)(dd);
                }
            }
        }
        // Load trees for the last iteration
        cur_tree   = param_tree.get_child_optional( "DT" );
        it_assert(cur_tree, "lut_ldpc::get_lut_tree_templates: Error reading decision node tree from file!");
        var_luts(max_iters-1).set_length(dv_act);
        for(int dd=0; dd<dv_act; dd++){
            ss << "var_deg_" << std::setfill('0') << std::setw(3) << var_deg(dd);
            boost::optional<std::string> cur_string = cur_tree->get_optional<std::string>( ss.str() );
            it_assert(cur_string, "lut_ldpc::get_lut_tree_templates(): Could not load tree string for variable node degree " << var_deg(dd) << " at iteration 0");
            var_luts(max_iters-1)(dd) = LUT_Tree(*cur_string, LUT_Tree::DECTREE);
            it_assert(var_luts(max_iters-1)(dd).get_num_leaves() == var_deg(dd)+1, "lut_ldpc::get_lut_tree_templates(): LUT tree does not match node degree!");
            var_luts(max_iters-1)(dd).set_resolution(Nq_Msg(0), Nq_Msg(1), Nq_Cha);
            ss.str(std::string());
            ss.clear();
        }
        
        if(!minLUT){
            // Initial Iteration
            cur_tree   = param_tree.get_child_optional( "chk_iter_000" );
            it_assert(cur_tree, "lut_ldpc::get_lut_tree_templates(): Error reading check node tree from file!");
            chk_luts.set_length(max_iters);
            chk_luts(0).set_length(dc_act);
            for(int dd=0; dd<dc_act; dd++){
                ss << "chk_deg_" << std::setfill('0') << std::setw(3) << chk_deg(dd);
                boost::optional<std::string> cur_string = cur_tree->get_optional<std::string>( ss.str() );
                it_assert(cur_string, "lut_ldpc::get_lut_tree_templates(): Could not load tree string for check node degree " << chk_deg(dd) << " at iteration 0");
                chk_luts(0)(dd) = LUT_Tree(*cur_string, LUT_Tree::CHKTREE);
                it_assert(chk_luts(0)(dd).get_num_leaves() == chk_deg(dd)-1, "lut_ldpc::get_lut_tree_templates(): LUT tree does not match node degree!");
                chk_luts(0)(dd).set_resolution(Nq_Msg(0), Nq_Msg(1));
                ss.str(std::string());
                ss.clear();
            }
            // Subsequent iterations
            for(int ii=1; ii< max_iters; ii++){
                chk_luts(ii).set_length(dc_act);
                ss << "chk_iter_" << std::setfill('0') << std::setw(3) << ii;
                cur_tree = param_tree.get_child_optional( ss.str() );
                ss.str(std::string());
                ss.clear();
                if(cur_tree){
                    for(int dd=0; dd<dc_act; dd++){
                        ss << "chk_deg_" << std::setfill('0') << std::setw(3) << chk_deg(dd);
                        boost::optional<std::string> cur_string  = cur_tree->get_optional<std::string>( ss.str() );
                        it_assert(cur_string, "lut_ldpc::get_lut_tree_templates(): Could not load tree string for check node degree " << chk_deg(dd) << " at iteration " << ii << ".");
                        chk_luts(ii)(dd) = LUT_Tree(*cur_string, LUT_Tree::CHKTREE);
                        it_assert(chk_luts(ii)(dd).get_num_leaves() == chk_deg(dd)-1, "lut_ldpc::get_lut_tree_templates(): LUT tree does not match node degree!");
                        chk_luts(ii)(dd).set_resolution(Nq_Msg(ii), Nq_Msg(ii+1), Nq_Cha);
                        ss.str(std::string());
                        ss.clear();
                    }
                }
                else{
                    for(int dd=0; dd<dv_act; dd++){
                        chk_luts(ii)(dd) = chk_luts(ii-1)(dd);
                    }
                }
            }
        }
        
    }
    else if( (tm=="auto_bin_balanced" || tm=="auto_bin_high" || tm=="root_only") && filename.length()==0){
        //Generate var tree structures
        var_luts.set_length(max_iters);
        for(int ii=0; ii< max_iters; ii++){
            var_luts(ii).set_length(dv_act);
            for(int dd=0; dd<dv_act; dd++){
                if(ii==max_iters-1){
                    var_luts(ii)(dd) = LUT_Tree(var_deg(dd)+1, LUT_Tree::DECTREE, tm);
                    var_luts(ii)(dd).set_resolution(Nq_Msg(ii), 2 , Nq_Cha);
                    it_info_debug("Generated Decision node tree for iteration " << ii << " and degree " << var_deg(dd)
                                  << ":" << var_luts(ii)(dd).gen_template_string());
                }
                else{
                    var_luts(ii)(dd) = LUT_Tree(var_deg(dd), LUT_Tree::VARTREE, tm);
                    var_luts(ii)(dd).set_resolution(Nq_Msg(ii), Nq_Msg(ii+1), Nq_Cha);
                    it_info_debug("Generated Var-LUT tree for iteration " << ii << " and degree " << var_deg(dd)
                                  << ":" << var_luts(ii)(dd).gen_template_string());
                }
            }
        }
        if(!minLUT){
            //Generate check tree structures
            chk_luts.set_length(max_iters);
            for(int ii=0; ii< max_iters; ii++){
                chk_luts(ii).set_length(dc_act);
                for(int dd=0; dd<dc_act; dd++){
                    chk_luts(ii)(dd) = LUT_Tree(chk_deg(dd)-1, LUT_Tree::CHKTREE, tm);
                    chk_luts(ii)(dd).set_resolution(Nq_Msg(ii), Nq_Msg(ii));
                    it_info_debug("Generated Chk-LUT tree for iteration " << ii << " and degree " << chk_deg(dd)
                                  << ":" << chk_luts(ii)(dd).gen_template_string());
                    
                }
            }
        }
        
    }
    else{
        it_error("Could Not parse tree_method " << tree_method);
    }
}


vec lut_ldpc::joint_level_irr_lut_design(const vec& degree_dist, const ivec& degrees, Array<LUT_Tree>& lut_trees, vec& P_row, double& Pe)
{
    int L = length(degree_dist);
    it_assert( L == lut_trees.length(), "joint_level_irr_lut_design(): Input dimension mismatch!");

    // Find the deepest level of the tree_array
    ivec levels(L);
    for(int ll=0; ll<L; ll++){
        levels(ll) = lut_trees(ll).get_height();
    //    lut_trees(ll).tikz_draw_tree("debtree" + to_str(ll) + ".tikz");
    }
    
    // Start processing the trees starting from their leaves
    int cur_level = max(levels) -1;
    
    vec avg_pmf;
    
    while(cur_level>=0){
        
        Array< deque<LUT_Tree_Node*> > level_nodes(L);
        LUT_Tree::tree_type_t t = lut_trees(0).get_type();
        for(int ll=0; ll<L; ll++){
            it_assert_debug(t == lut_trees(ll).get_type() , "joint_level_irr_lut_design(): input trees have inconsistent type");
            if(levels(ll) > cur_level){
                level_nodes(ll) = lut_trees(ll).get_level_nodes(cur_level);
                // Remove nodes of unwanted type and get node weights
                deque<LUT_Tree_Node*>::iterator ii;
                for ( ii = level_nodes(ll).begin(); ii != level_nodes(ll).end() ; )                {
                    if ( ((*ii)->type == LUT_Tree_Node::IM) || ((*ii)->type == LUT_Tree_Node::ROOT) )
                        ++ii; //  keep element and proceed
                    else{
                        ii = level_nodes(ll).erase(ii); // erase returns the next iterator
                    }
                }
            }
        }
        // Assemble nodes and weights of current level and perform the LUT design
        avg_pmf = level_lut_tree_update(level_nodes, degree_dist, t );
        cur_level--;
    }
    if(length(P_row)>0){
        Pe = 0;
        for(int dd=0; dd<L; dd++){
            vec p_tmp = lut_trees(dd).update(true);
            P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
            Pe += degree_dist(dd)*P_row(dd);
        }
    }
    
    return avg_pmf;
}

vec lut_ldpc::joint_root_irr_lut_design(const vec& degree_dist, const ivec& degrees, Array<LUT_Tree>& lut_trees, vec& P_row, double& Pe)
{
    int L = length(degree_dist);
    it_assert( L == lut_trees.length(), "joint_root_irr_lut_design(): Input dimension mismatch!");
    
    for(int dd=0; dd<L; dd++){
        (void) lut_trees(dd).update();
    }
    
    vec avg_pmf;

    
    Array< deque<LUT_Tree_Node*> > root_nodes(L);
    LUT_Tree::tree_type_t t = lut_trees(0).get_type();
    for(int ll=0; ll<L; ll++){
        it_assert_debug(t == lut_trees(ll).get_type() , "joint_level_irr_lut_design(): input trees have inconsistent type");
            root_nodes(ll) = lut_trees(ll).get_level_nodes(0);
        }
    
    // Assemble nodes and weights of current level and perform the LUT design
    avg_pmf = level_lut_tree_update(root_nodes, degree_dist, t );

    if(length(P_row)>0){
        Pe = 0;
        for(int dd=0; dd<L; dd++){
            vec p_tmp = lut_trees(dd).update(true);
            P_row(dd) = sum(p_tmp.left(p_tmp.size()/2));
            Pe += degree_dist(dd)*P_row(dd);
        }
    }
    
    return avg_pmf;
}

vec lut_ldpc::level_lut_tree_update(Array< deque<LUT_Tree_Node*> >& tree_nodes,  const vec& degree_dist, LUT_Tree::tree_type_t t)
{
    int L = tree_nodes.size();
    it_assert(L == length(degree_dist), "level_lut_tree_update(): Dimension mismatch");
    // Get the node weights and product pmfs
    Array<vec> node_weights(L);
    Array<Array<vec>> pmf_prod(L);
    Array<ivec> pmf_prod_len(L);
    int M_tot = 0; //overall number of input labels for quantizer design
    int Num_outlabels = -1;
    
    for(int ll=0; ll<L; ll++){
        int J = (int) tree_nodes(ll).size();
        pmf_prod(ll).set_size(J);
        pmf_prod_len(ll).set_size(J);
        node_weights(ll)= vec(J);
        for(int jj=0; jj<J; jj++){
            if (Num_outlabels == -1)  Num_outlabels =  tree_nodes(ll)[jj]->K;
            it_assert_debug(Num_outlabels == tree_nodes(ll)[jj]->K, "level_lut_tree_update(): Output resolution mismatch");
            node_weights(ll)(jj) = tree_nodes(ll)[jj]->get_num_leaves();
            pmf_prod(ll)(jj) = tree_nodes(ll)[jj]->get_input_product_pmf(t);
            pmf_prod_len(ll)(jj) = pmf_prod(ll)(jj).length();
            M_tot += pmf_prod_len(ll)(jj);
        }
        node_weights(ll) = node_weights(ll)/sum(node_weights(ll));
    }
    
    
    vec pmf_prod_overall = -1e9*ones(M_tot);
    int I = 0;
    for(int ll=0; ll<L; ll++){
        int J = (int) tree_nodes(ll).size();
        for(int jj=0; jj<J; jj++){
            int M =  pmf_prod_len(ll)(jj);
            for(int mm=0; mm<M/2; mm++){
                pmf_prod_overall(I + mm) = node_weights(ll)(jj)*degree_dist(ll)*pmf_prod(ll)(jj).get(mm);
                pmf_prod_overall(M_tot-1 -I - mm) = node_weights(ll)(jj)*degree_dist(ll)*pmf_prod(ll)(jj).get(M-1-mm);
            }
            I += M/2;
        }
    }

//     The pmf does not necessarily sum to 1 as not all degrees might be active for high layers
//    cout << "pmf_prod_overall = " << pmf_prod_overall << endl << sum(pmf_prod_overall) << endl << endl;
    pmf_prod_overall = pmf_prod_overall/sum(pmf_prod_overall);
    
    // Joint Quantizer design for root nodes
    ivec Q_out_overall;
    vec p_out;
    // Eliminate entries with mass 0
    bvec idx_nz = (.5*(pmf_prod_overall + fliplr(pmf_prod_overall)) != 0);
    ivec Q_out_overall_nz;
    (void) quant_mi_sym(p_out, Q_out_overall_nz, pmf_prod_overall.get(idx_nz), Num_outlabels);
    // For the entries with mass zero we set the outputs symmetrically and assign the least confident llr magnitudes
    Q_out_overall = concat(ones_i(M_tot/2)*(Num_outlabels/2-1), ones_i(M_tot/2)*(Num_outlabels/2) );
    int mm_nz=0;
    for(int mm=0; mm<M_tot; mm++){
        if(idx_nz(mm)){
            Q_out_overall(mm) = Q_out_overall_nz(mm_nz);
            mm_nz++;
        }
    }
    
    Array<Array<ivec>> luts(L);
    I = 0;
    for(int ll=0; ll<L; ll++){
        int J = (int) tree_nodes(ll).size();
        luts(ll).set_size(J);
        for(int jj=0; jj<J; jj++){
            int M =  pmf_prod_len(ll)(jj);
            luts(ll)(jj).set_size(M/2);
            for(int mm=0; mm<M/2; mm++){
                luts(ll)(jj)(mm) = Q_out_overall(I+mm);
            }
            I += M/2;
            // Update the tree node quantizer map and output distribution
            tree_nodes(ll)[jj]->Q = luts(ll)(jj);
            tree_nodes(ll)[jj]->p = zeros(Num_outlabels);
            for(int mm=0; mm<M; mm++){
                if(mm<M/2)
                    tree_nodes(ll)[jj]->p(luts(ll)(jj)(mm)) += pmf_prod(ll)(jj)(mm);
                else
                    tree_nodes(ll)[jj]->p(Num_outlabels-1-luts(ll)(jj)(M-1-mm)) += pmf_prod(ll)(jj)(mm);
            }
        }
    }
    return p_out;
}





double lut_ldpc::get_lam2stable_qbp(double sig, vec rho, int Nq_Cha, double LLR_max, int Nq_fine){
    double delta = 2*LLR_max/Nq_fine;
    vec pmf_channel_fine = get_gaussian_pmf(2/sqr(sig), 2/sig, Nq_fine, delta);
    
    // Do not consider degree one
    rho.del(0);
    ivec Q_out;
    vec pmf_cha;
    (void) quant_mi_sym(pmf_cha,     Q_out, pmf_channel_fine, Nq_Cha, true);
    
    
    double e_to_r = 1 / sum(sqrt(elem_mult(pmf_cha, fliplr(pmf_cha))));
    double len_rho = length(rho);
    double rho_dev_eval_one = sum(elem_mult(rho, linspace(1, len_rho, len_rho)));
    return e_to_r / rho_dev_eval_one;
}

double lut_ldpc::get_lam2stable_cbp(double sig, vec rho){
    rho.del(0);
    int len_rho = rho.length();
    double rho_derivative_eval_0 =  sum(elem_mult(rho, linspace(1, len_rho, len_rho))); 
    return std::exp(1/(2*itpp::sqr(sig)))/rho_derivative_eval_0;
}

double lut_ldpc::get_lam2stable_qbp_iterative(double sig, vec rho, int Nq_Cha, double LLR_max, int Nq_fine){
    
    int Nbit = 13;
    int N = pow2i(Nbit-1);
    int Imax = 1e5;
    double cauchy_interval = 1e-9;
    
    // === Prepare channel pmf
    double delta = LLR_max/N;
    vec pmf_channel_fine = get_gaussian_pmf(2/sqr(sig), 2/sig, 2*N+2, delta);
    ivec Q_out;
    vec pmf_cha;
        // Get quantized channel pmf at high resolution
    (void) quant_mi_sym(pmf_cha,     Q_out, pmf_channel_fine, Nq_Cha, true);
    int ll = 0;
    vec pmf_cha_sparse = zeros(2*N+2);
    for(int nn=0; nn<2*N+1; nn++){
        double L = std::log(pmf_cha(ll)) - std::log(pmf_cha(Nq_Cha-1-ll));
        int nn_signed = nn - N;
        if(L > nn_signed*delta && L <=(nn_signed+1)*delta ){
            pmf_cha_sparse(nn) = pmf_cha(ll);
            ll++;
            if(ll>=Nq_Cha) break;
        }
    }
    
    int Nfft = pow2i(1+std::ceil(std::log2(2*N+1)));
    
    vec pmf_in = pmf_cha_sparse;
    vec pmf_out = pmf_cha_sparse;
    
//    vec pmf_in = pmf_channel_fine;
//    vec pmf_out = pmf_channel_fine;
    
    vec var_conv_weight = exp(-.5* linspace_fixed_step(-N, N)*delta);
    
    
    // Prepare CN distribution
    rho.del(0);
    double len_rho = length(rho);
    double rho_dev_eval_one = sum(elem_mult(rho, linspace(1, len_rho, len_rho)));

    double e_to_r;
    double e_to_r_old = std::numeric_limits<double>::min();
    
    ofstream myfile("lam2_trace_qbp4bit_iter_rich.txt");
    for(int ii=2; ii<Imax; ii++){
        it_assert(Nfft > 2*N+2, "SMLDPC_DE_PMF::var_update_convolve(): FFT size too small");
        // Scale and make even pmfs
//        vec a_hat =  elem_mult(var_conv_weight, pmf_in.left(2*N+1));
//        vec b_hat =  elem_mult(var_conv_weight, pmf_out.left(2*N+1));
        vec a_hat =  pmf_in.left(2*N+1);
        vec b_hat =  pmf_out.left(2*N+1);
        // convolve in frequency domain
        vec pmf_tmp = ifft_real( elem_mult( fft_real(a_hat, Nfft), fft_real(b_hat,Nfft) ) );
        
        // Cut out part correspondin to input pmf's support, neglecting the probability mass in [-2*N,-(N+1)] We obtain mass in [N+1, 2N] later
        pmf_out = pmf_tmp.get(N, 3*N);
        // Undo scaling
//        pmf_out = elem_div(pmf_tmp, var_conv_weight);
        // Compute mass at -delta*n and inf
        pmf_out(0) += sum(pmf_tmp.get(0, N-1));
        pmf_out = concat(pmf_out, 1 - sum(pmf_out));
        
  
        double Pe = sum(pmf_out.get(0, N-1))+ .5*pmf_out(N);
        if(Pe == 0) break;
        e_to_r = std::exp( - std::log(Pe)/ii);
        if( fabs(e_to_r_old - e_to_r) < cauchy_interval){
            break;
        }
        e_to_r_old = e_to_r;
        
        myfile << e_to_r / rho_dev_eval_one << endl;
    }
    myfile.close();
    return e_to_r / rho_dev_eval_one;
}

double lut_ldpc::get_lam2stable_lut(double sig, vec rho, int Nq_Cha, int Nq_Msg, double LLR_max, int Nq_fine){
    double delta = 2*LLR_max/Nq_fine;
    vec pmf_channel_fine = get_gaussian_pmf(2/sqr(sig), 2/sig, Nq_fine, delta);
    
    // Do not consider degree one
    rho.del(0);
    ivec Q_out;
    vec pmf_cha;
    vec pmf_con;
    (void) quant_mi_sym(pmf_cha,     Q_out, pmf_channel_fine, Nq_Cha, true);
    (void) quant_mi_sym(pmf_con,     Q_out, pmf_cha, Nq_Msg, true);
    
    double e_to_r;
    double e_to_r_old = std::numeric_limits<double>::min();
    int Nmax = 1e5;
    double cauchy_interval = 1e-6;
    
    for(int nn=0; nn<Nmax; nn++){
                
        Array <vec> p_in(2);
        p_in(0) = pmf_con;
        p_in(1) = pmf_cha;
        vec p_Msg_prod = get_var_product_pmf(p_in);
        // Eliminate entries with mass 0
        bvec idx_nz = (.5*(p_Msg_prod + fliplr(p_Msg_prod)) != 0);
        ivec Q_out_nz;
        (void) quant_mi_sym(pmf_con, Q_out_nz, p_Msg_prod.get(idx_nz), Nq_Msg);
      //  e_to_r = std::exp( - std::log(sum(pmf_con.left(Nq_Msg/2)))/nn);
        e_to_r = std::pow(sum(pmf_con.left(Nq_Msg/2)), -1.0/nn);

        if( fabs(e_to_r_old - e_to_r) < cauchy_interval){
            break;
        }
        e_to_r_old = e_to_r;
    }
    
    double len_rho = length(rho);
    double rho_dev_eval_one = sum(elem_mult(rho, linspace(1, len_rho, len_rho)));
    return e_to_r / rho_dev_eval_one;
}

std::ostream& lut_ldpc::operator<<(std::ostream &os, const LDPC_Ensemble &ens){
    TextTable l( '-', '|', '+' );
    l.add( "VN degrees" );
    for(int ii=0; ii<ens.dv_act; ii++){
        std::ostringstream strs;
        strs << ens.degree_lam(ii);
        l.add( strs.str() );
    }
    l.endOfRow();
    l.add( "VN edge pmf" );
    for(int ii=0; ii<ens.dv_act; ii++){
        std::ostringstream strs;
        strs << ens.lam(ii);
        l.add( strs.str() );
    }
    l.endOfRow();
    TextTable r( '-', '|', '+' );
    r.add( "CN degrees" );
    for(int ii=0; ii<ens.dc_act; ii++){
        std::ostringstream strs;
        strs << ens.degree_rho(ii);
        r.add( strs.str() );
    }
    r.endOfRow();
    r.add( "CN edge pmf" );
    for(int ii=0; ii<ens.dc_act; ii++){
        std::ostringstream strs;
        strs << ens.rho(ii);
        r.add( strs.str() );
    }
    r.endOfRow();
    
    os << l << r;
    return os;
}


