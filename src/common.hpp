/*!
 * \file common.hpp
 * \brief Collection of function that are accesses throughout lut_ldpc
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

//! Nonlinear quantization of \c x. The output is equal to i if \c boundaries(i-1) < x < boundaries(i)
int quant_nonlin(double x, const vec& boundaries);
//! Vector version of quant_nonlin
ivec quant_nonlin(const vec& x, const vec& boundaries);
    
//! Linear quantization of x
int quant_lin(double x, double delta, int N);

/*! \brief Get quantized Gaussian pmfs with N quantization intervals (2 Overload and
 N -2) inner regions. For N odd, this is quantization with 0, with N even there is no zero.
 */
vec get_gaussian_pmf(double mu, double sig, int N, double delta);

//! Convert rate to noise threshold
double rate_to_shannon_thr(double R);
//! Get maximum rate for noise theshold
double shannon_thr_to_rate(double sig);

//! x*log2(y) with checking for numeric stability
inline double x_log2_y(double x, double y);

//! Flip contents of a vector, i.e for a length N vector,  x_out(ii) = x_in(N-1 -ii)
template <class Num_T> Vec<Num_T> fliplr(const Vec<Num_T>& x);

//! Eliminate duplicate entries from a vector
template <class Num_T> Vec<Num_T> unique(const Vec<Num_T>& x);

//! Kronecker product of two vectors
template <class Num_T> Vec<Num_T> kron(const Vec<Num_T>& x, const Vec<Num_T>& y);

//! Map sign based index to magnitude based index
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
