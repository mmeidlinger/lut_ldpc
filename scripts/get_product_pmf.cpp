#include <iostream>
#include <math.h>
#include "mex.h"
#include <itpp/itbase.h>
#include <itpp/itmex.h>


#include "../src/LDPC_DE.hpp"

// compiled  with
// mex  -I/usr/local/include -I"/Users/mmeidlin/Work/epfl-tuwien-ldpc/cpp/osx/include" -L"/Users/mmeidlin/Work/epfl-tuwien-ldpc/cpp/osx/lib" get_product_pmf.cpp ../src/LDPC_DE.cpp  -litpp_static_debug -llapack -lfftw3 -lblas -lboost_system -lboost_filesystem

using namespace itpp;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    // Check the number of inputs and output arguments
    if(nlhs!=1)  mexErrMsgTxt("Wrong number of output variables!");
    if(nrhs > 4)  mexErrMsgTxt("Wrong number of input variables!: Usage: get_chk_prod_pmf(dv,p, mode, Nq_optional)");
    
    int d = mxArray2int(prhs[0]);
    vec p_in= mxArray2vec(prhs[1]);
    int Nq;    
    vec p_out;

    // Get Product pmf
    if(nrhs == 3 ){ // Check node
        Nq = mxArray2int(prhs[2]);
        Array<vec> p_array(d-1);
        for(int ii=0; ii< d-1; ii++){
            p_array(ii)=p_in;   
        }
        p_out = get_chk_product_pmf(p_array);
    }
    else if(nrhs == 4){ // Variable node
        vec p_cha = mxArray2vec(prhs[2]);
        Nq = mxArray2int(prhs[3]);
        Array<vec> p_array(d);
        for(int ii=0; ii< d-1; ii++){
            p_array(ii)=p_in;   
        }
        p_array(d-1) = p_cha;
        p_out = get_var_product_pmf(p_array);
    }
    else{
        mexErrMsgTxt("nrhs must be 3 (check nodes) or 4 (variable nodes)");
    }
    
    // Eliminate Duplicate LLRs
    ivec idx1, idx2;
    p_out = sym_llr_sort_unique(p_out, idx1, idx2, 1e-10);
    
    // Quantize if necessary
    if(Nq > 0){
        vec p_tmp = p_out; 
        ivec Q;
        (void) quant_mi_sym(p_out, Q, p_tmp, Nq, true);
    }
    
    plhs[0] = mxCreateDoubleMatrix(1,p_out.size(), mxREAL);
    vec2mxArray(p_out, plhs[0]);
    return;
}
