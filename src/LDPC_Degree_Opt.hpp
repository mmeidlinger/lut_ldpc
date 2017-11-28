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



#ifndef LDPC_Degree_Opt_hpp
#define LDPC_Degree_Opt_hpp

#include "LDPC_DE.hpp"
#include <glpk.h>

//! Git version of program
extern const char *gitversion;

//! @namespace itpp IT++ Namespace
namespace itpp{
// ----------------------------------------------------------------------
// LDPC_Degree_Opt_LP
// ----------------------------------------------------------------------

//! Class for degree distribution optimization
class LDPC_Degree_Opt_LP
{
    //https://math.stackexchange.com/questions/432003/converting-absolute-value-program-into-linear-program
public:
    //! Constrict LP object
    LDPC_Degree_Opt_LP(LDPC_DE* de, double rate_target,
                               bool var_opt,
                               bool chk_opt,
                               double delta,
                               double scale_down, 
                               bool lam2constraint,
                               int maxiter_lp):
        den_evol(de), rate_target(rate_target), var_opt(var_opt),
        chk_opt(chk_opt), delta(delta), scale_down(scale_down), lam2constraint(lam2constraint),
        maxiter_lp(maxiter_lp), lam2stable(1.0) {};

    //! Run a degree optimization using an iterative, linear program, local hill clambing method
    int run(LDPC_Ensemble& ens);
    //! Solves the linear program needed for the hill climbing method
    void opt_lp(LDPC_Ensemble& ens, const mat& P, const vec& p, int mode);
    
    
private:
    
    
    //! Density evolution
    LDPC_DE* den_evol;
    //! Rate target
    double rate_target;
    //! True if variable node degree distribution is subject to optimization
    bool var_opt;
    //! True if check node degree distribution is subject to optimization
    bool chk_opt;
    //! LP closeness constraint value/stepsize
    double delta;
    //! If an updated ensemble does not converge, we scale down the channel standard deviation by this factor    
    double scale_down;
    
    
    //! Wether a constraint on variable node degree 2 is placed on LP optimization
    bool lam2constraint;

    //! Maximum number of LP iterations
    int maxiter_lp;
        
    //! Maximum stable variable node edge degree 2
    double lam2stable;

    //! Wether variable or check node distributions are being optimized
    enum{
        VAR_OPT,
        CHK_OPT
    };
    
};

}
#endif
