#define LLR_MAX 25.0
#define NQ_FINE 5000

#include <iostream>
#include <math.h>
#include "mex.h"
#include <itpp/itbase.h>
#include <itpp/itmex.h>


#include "../src/LDPC_DE.hpp"

// compiled  with
// mex  -I/usr/local/include -I"/Users/mmeidlin/Work/epfl-tuwien-ldpc/cpp/osx/include" -L"/Users/mmeidlin/Work/epfl-tuwien-ldpc/cpp/osx/lib" get_quant_biawgn_pmf.cpp ../src/LDPC_DE.cpp  -litpp_static_debug -llapack -lfftw3 -lblas -lboost_system -lboost_filesystem

using namespace itpp;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    // Check the number of inputs and output arguments
    if(nlhs!=1)  mexErrMsgTxt("Wrong number of output variables!");
    if(nrhs!=2)  mexErrMsgTxt("Wrong number of input variables!: Usage: get_quant_biawgn_pmf(sig2, Nq)");
    
    double sig2 = mxArray2double(prhs[0]);
    int Nq = mxArray2double(prhs[1]);
    double delta = 2*LLR_MAX/NQ_FINE;
    double sig = std::sqrt(sig2);
    //double delta = 2*(2/sqr(sig) + 4*2/sig)/Nq_fine;
    vec pmf_channel_fine = get_gaussian_pmf(2/sqr(sig), 2/sig, NQ_FINE, delta);
    
    ivec Q_out;
    vec pmf_cha;
    (void) quant_mi_sym(pmf_cha,     Q_out, pmf_channel_fine, Nq, true);
    plhs[0] = mxCreateDoubleMatrix(1,pmf_cha.size(), mxREAL);
    vec2mxArray(pmf_cha, plhs[0]);
    return;
}
