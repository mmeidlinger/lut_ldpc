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

#include "LDPC_Degree_Opt.hpp"
using namespace itpp;
using namespace std;


int LDPC_Degree_Opt_LP::run(LDPC_Ensemble& ens){
    
    int ii;
    double thr;
    mat P;
    vec p;
    
    den_evol->set_ensemble(ens);
    double thr_max = rate_to_shannon_thr(ens.get_rate());
    double thr_min = thr_max/1e4;
    den_evol->set_bisec_window(thr_min, thr_max);
    
    
    (void) den_evol->bisec_search(thr);
    cout << "Initial threshold = " << thr << ", initial Rate = " << ens.get_rate() << endl; 
    
    lam2stable = std::max(den_evol->get_lam2stable(thr), ens.get_lam_of_degree(2)); 
    
    for(ii=0; ii< maxiter_lp; ii++){
        
        // optimize variable node degree
        if(var_opt){
            int iters  = den_evol->evolve(thr, 1, 0,  P, p);
            while (iters <= 0) {
                thr *= scale_down;
                cout << "Pre VN optimization: Updated ensemble" << endl << ens << "does not converge,"
                         " scaling down channel paramenter sigma to " << thr << endl;
                iters = den_evol->evolve(thr, 1, 0,  P, p);
            }
            p =  P*ens.sget_lam();
            
            
            opt_lp(ens, P, p, VAR_OPT);
            den_evol->set_ensemble(ens);
            //cout << "   Optimized var: Rate = " << ens.get_rate() << " thr = " << thr << endl;
        }
        
        // optimize check node degree
        if(chk_opt){
            int iters = den_evol->evolve(thr, 0, 1,  P, p);
            while (iters <= 0) {
                thr *= scale_down;
                cout << "Pre CN optimization: Updated ensemble" << endl << ens << "does not converge,"
                             " scaling down channel paramenter sigma to " << thr << endl;
                iters = den_evol->evolve(thr, 0, 1,  P, p);
            }
            p =  P*ens.sget_rho();
            opt_lp(ens, P, p, CHK_OPT);
            den_evol->set_ensemble(ens);
            //cout << "   Optimized  chk: Rate = " << ens.get_rate() << " thr = " << thr << endl;
        }
        
        lam2stable = std::max(den_evol->get_lam2stable(thr), ens.get_lam_of_degree(2));
        
        cout << "== LP-ITERATION " << ii << endl;
        cout << "Rate = " << ens.get_rate() << " thr = " << thr << endl;
        cout << "Var Edge distribution = " << endl << ens.sget_degree_lam() << endl << ens.sget_lam() << endl;
        cout << "Chk Edge distribution = " << endl << ens.sget_degree_rho() << endl << ens.sget_rho() << endl;
        cout << "Maximum Stable VN degree = " << lam2stable << endl << endl;
        // check for convergence
        if(ens.get_rate() > rate_target ){
            break;
        }
    }
    return ii;
}

void LDPC_Degree_Opt_LP::opt_lp(itpp::LDPC_Ensemble &ens, const mat &P, const vec &p, int mode){
    it_assert( P.rows() == p.length(), "SMLDPC_Degree_Opt_LP::opt_lp(): Input Dimension Mismatch");
    
    int iters = P.rows(); 
    int dim;
    ivec deg;
    double Pe0 = 1.0;
    
    switch (mode) {
        case VAR_OPT:
            dim = ens.get_dv_act();
            deg = ens.sget_degree_lam();
            break;
        case CHK_OPT:
            dim = ens.get_dc_act();
            deg = ens.sget_degree_rho();
            break;
        default:
            it_error("SMLDPC_Degree_Opt_LP::var_opt_sym(): Mode not defined");
            break;
            
    }
    it_assert(dim == P.cols(), "SMLDPC_Degree_Opt_LP::sym_chk_opt_lp(): Input Dimension Mismatch");
    
    // Build linear program and allocate memory-
    glp_prob *lp;
    lp = glp_create_prob();
    
    switch (mode) {
        case VAR_OPT:
            glp_set_obj_dir(lp, GLP_MAX);
            break;
        case CHK_OPT:
            glp_set_obj_dir(lp, GLP_MIN);
            break;
        default:
            it_error("SMLDPC_Degree_Opt_LP::var_opt_sym(): Mode not defined");
            break;
    }

    
    
    
    int num_matrix_el = dim* (3*iters+1);
    int    *ia = (int*)    malloc(sizeof(int)    * (num_matrix_el+1));
    int    *ja = (int*)    malloc(sizeof(int)    * (num_matrix_el+1));
    double *ar = (double*) malloc(sizeof(double) * (num_matrix_el+1));
    
    if( (ia && ja && ar) == false){
       it_error("SMLDPC_Degree_Opt_LP::var_opt_sym(): Could not allocate memory");
    }
    // Set columns (= variables to optimize over)
    glp_add_cols(lp, dim);
    for(int ii=0; ii<dim; ii++){
        std::ostringstream col_name;
        col_name << "Edge degree probability for degree " << deg(ii) << " edges";
        glp_set_col_name(lp, ii+1, col_name.str().c_str());
        // Set both upper and lower bounds for degree 2, lower bound 0 for all other degrees
        if(ii==0 && mode==VAR_OPT && lam2constraint) glp_set_col_bnds(lp, ii+1, GLP_DB, 0.0, lam2stable);
        else glp_set_col_bnds(lp, ii+1, GLP_LO, 0.0, 0.0); // lower bound 0 (second zero is ignored)
        glp_set_obj_coef(lp, ii+1, 1.0/ (double)deg(ii));
    }
    
    // Set rows (= constraints)
    glp_add_rows(lp, 1+3*iters);
    int mm=1;
    
    glp_set_row_name(lp, 1, "sum one constraint");
    glp_set_row_bnds(lp, 1, GLP_FX, 1.0, 1.0);
    for(int ii=0; ii<dim; ii++){
        ia[mm] = 1;
        ja[mm] = ii+1;
        ar[mm] = 1;
        mm++;
    }
    
    
    glp_set_row_name(lp, 2, "closeness constraint, pos, iteration=0");
    glp_set_row_bnds(lp, 2, GLP_UP, 0.0,  p(0)+fmax(0, delta*(Pe0-p(0))));
    for(int ii=0; ii<dim; ii++){
        ia[mm] = 2;
        ja[mm] = ii+1;
        ar[mm] = P(0,ii);
        mm++;
    }
    
    glp_set_row_name(lp, 3, "closeness constraint, neg, iteration=0");
    glp_set_row_bnds(lp, 3, GLP_UP, 0.0, -p(0)+fmax(0, delta*(Pe0-p(0))));
    for(int ii=0; ii<dim; ii++){
        ia[mm] = 3;
        ja[mm] = ii+1;
        ar[mm] = -P(0,ii);
        mm++;
    }
    
    glp_set_row_name(lp, 4, "improvement constraint,    iteration=0");
    glp_set_row_bnds(lp, 4, GLP_UP, 0.0, Pe0);
    for(int ii=0; ii<dim; ii++){
        ia[mm] = 4;
        ja[mm] = ii+1;
        ar[mm] = P(0,ii);
        mm++;
    }
    
    for(int ll=1; ll<iters; ll++){
        
        std::ostringstream row_name;
        
        row_name << "closeness constraint, pos, iteration=" << ll;
        glp_set_row_name(lp, 3*ll+2, row_name.str().c_str());
        glp_set_row_bnds(lp, 3*ll+2, GLP_UP, 0.0,  p(ll)+fmax(0, delta*(p(ll-1)-p(ll))));
        for(int ii=0; ii<dim; ii++){
            ia[mm] = 3*ll+2;
            ja[mm] = ii+1;
            ar[mm] = P(ll,ii);
            mm++;
        }
        row_name.str("");
        row_name.clear();
        
        row_name << "closeness constraint, neg, iteration=" << ll;
        glp_set_row_name(lp, 3*ll+3, row_name.str().c_str());
        glp_set_row_bnds(lp, 3*ll+3, GLP_UP, 0.0, -p(ll)+fmax(0, delta*(p(ll-1)-p(ll))));
        for(int ii=0; ii<dim; ii++){
            ia[mm] = 3*ll+3;
            ja[mm] = ii+1;
            ar[mm] = -P(ll,ii);
            mm++;
        }
        row_name.str("");
        row_name.clear();
        
        row_name << "improvement constraint,    iteration=" << ll;
        glp_set_row_name(lp, 3*ll+4, row_name.str().c_str());
        glp_set_row_bnds(lp, 3*ll+4, GLP_UP, 0.0, p(ll-1));
        for(int ii=0; ii<dim; ii++){
            ia[mm] = 3*ll+4;
            ja[mm] = ii+1;
            ar[mm] = P(ll,ii);
            mm++;
        }
    }
    
    it_assert(mm == num_matrix_el+1, "SMLDPC_Degree_Opt_LP::var_opt_sym(): Internal Error");
    
    // setup and solve lp
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
    
    glp_load_matrix(lp, num_matrix_el, ia, ja, ar);
    int solver_return_value = glp_simplex(lp, &parm);
    //glp_print_sol(lp, "glpk_symchk_sol.txt");
    
    // write solution to ensemble
    vec solution_vec(dim);
    for(int ii=0; ii<dim; ii++){
        solution_vec(ii) = glp_get_col_prim(lp, ii+1);
    }
    // Warning: Return value 0 does not imply that the solution is feasible!
    if(solver_return_value == 0){
        switch (mode) {
            case VAR_OPT:
                ens.sset_lam(solution_vec);
                break;
            case CHK_OPT:
                ens.sset_rho(solution_vec);
                break;
            default:
                it_error("SMLDPC_Degree_Opt_LP::var_opt_sym(): Mode not defined");
                break;
        }
    }
    // clean up
    glp_delete_prob(lp);
    if(ia) free(ia), ia = NULL;
    if(ja) free(ja), ja = NULL;
    if(ar) free(ar), ar = NULL;
    
    return;
}
