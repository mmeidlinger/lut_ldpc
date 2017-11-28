//
//  common.hpp
//  LDPC_LUT
//
//  Created by Michael Meidlinger on 28.11.17.
//

#ifndef common_hpp
#define common_hpp

//! LLRs are considered identical within this tolerance
#define UNIQE_LLR_DELTA 0.0

#include <stdio.h>
#include <itpp/itbase.h>

using namespace itpp;


namespace lut_ldpc{

//! Get the product distribution for a variable node with input probabilities p_in
vec get_var_product_pmf(const Array<vec>& p_in);

//! Get the product distribution for a check node with input probabilities p_in
vec get_chk_product_pmf(const Array<vec>& p_in);
    
//! Convert SNR (=Eb/N0 in dB) to corresponding AWGN channel noise standard deviation
inline double snr2sig(double rate, double snr);
//! Convert AWGN channel noise standard deviation to SNR (=Eb/N0 in dB)
inline double sig2snr(double rate, double sig);

//! Convert AWGN channel noise standard deviation to SNR (=Eb/N0 in dB)
vec sig2snr(double rate, const vec& sig);

//! Convert SNR (=Eb/N0 in dB) to corresponding AWGN channel noise standard deviation
vec snr2sig(double rate, const vec& snr);

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
    
template <class Num_T> ivec sort_index_sym(const Vec<Num_T>& x);
template<class Num_T> ivec get_complement_idx(const Vec<Num_T>& a, const Vec<Num_T>& b,  Num_T c);

int signed_to_unsigned_idx(int idx, const ivec& inres);
    
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
 * @returns     Returns a symmetric pmf with unique LLRs and the corresponding indices
 * @param[in]   p_in        An arbitrary, unsorted, conditional pmf.
 * @param[out]  idx_in      sorting index of input pmf according to LLR
 * @param[out]  idx_sorted  index mapping input to unique pmf
 */
vec sym_llr_sort_unique(const vec& p_in, ivec& idx_in, ivec& idx_sorted, double llr_delta = UNIQE_LLR_DELTA);

//! Calculate the mutual information between X and Y, where p(y|x)=p(-y|-x) is given by p_in and X is binary and uniform
double get_mi_bcpmf_sym(const vec& p);
    
}// namespace lut_ldpc

#endif /* common_hpp */
