#include <itpp/itcomm.h>
#include <itpp/itmex.h>
#include <itpp/itbase.h>

// compiled  with
// mex  -I~/Work/epfl-tuwien-ldpc/cpp/include matlab_matrix_to_alist.cpp -L~/Work/epfl-tuwien-ldpc/cpp/lib -litpp_static -llapack -lfftw3 -lblas

using namespace itpp;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    // Check the number of inputs and output arguments
    if(nlhs!=0)  mexErrMsgTxt("Wrong number of output variables!");
    if(nrhs!=5)  mexErrMsgTxt("Wrong number of input variables!: Usage: sparse_matlab_matrix_to_alist(num_rows, num_cols, row_idx_vec, col_idx_vec, outname.alist )");
    
    int num_rows = mxArray2int(prhs[0]);
    int num_cols = mxArray2int(prhs[1]);
    ivec row_idx = mxArray2ivec(prhs[2]);
    ivec col_idx = mxArray2ivec(prhs[3]);
    std::string filename = mxArray2string(prhs[4]);
    
    it_assert(num_rows>0 && num_cols>0, "Number of rows and number of columns must be positive integers!");
    it_assert(row_idx.length() == col_idx.length(), "Row index and column index vector must have same length!");
    
    int num_edges = row_idx.length();
    GF2mat_sparse H(num_rows, num_cols);
    
    for(int ee=0; ee<num_edges; ee++){
        H.set(row_idx(ee),col_idx(ee), 1);
    }
    
    GF2mat_sparse_alist Halist;
    Halist.from_sparse(H);
    Halist.write(filename);
    
    return;
}
