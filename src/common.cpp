/*!
 * \file common.cpp
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

#include "common.hpp"
vec lut_ldpc::get_var_product_pmf(const Array<vec>& p_in){
    // Build product distribution
    int num_inputs = p_in.size();
    vec p_Msg_prod = p_in(num_inputs-1);
    
    for(int ii=num_inputs-2; ii>=0; ii--){
        p_Msg_prod = kron(p_Msg_prod, p_in(ii));
    }
    return p_Msg_prod;
}

vec lut_ldpc::get_chk_product_pmf(const Array<vec>& p_in){
    // Build product distribution
    int num_inputs = p_in.size();
    ivec res_inputs(num_inputs);
    for(int jj=0; jj<num_inputs; jj++){
        res_inputs(jj) = p_in(jj).length();
    }
    
    
    vec p_Msg_prod0 = p_in(num_inputs-1);
    vec p_Msg_prod1 = fliplr(p_in(num_inputs-1));
    
    // In this loop, in order to respect the label order, we have p_Msg_prod0  != fliplr(p_Msg_prod1) in general. Symmetry is restored by mapping signed to unsigned labels
    for(int ii=num_inputs-2; ii>=0; ii--){
        vec p_Msg_prod0_new = .5* (kron(p_Msg_prod0, p_in(ii))  + kron(p_Msg_prod1, fliplr(p_in(ii))));
        vec p_Msg_prod1_new = .5* (kron(p_Msg_prod1, p_in(ii))  + kron(p_Msg_prod0, fliplr(p_in(ii))));
        p_Msg_prod0 = p_Msg_prod0_new;
        p_Msg_prod1 = p_Msg_prod1_new;
    }
    
    vec p_Msg_prod_combined = zeros( 2*prod(res_inputs/2));
    //re-order and combine labels
    for (int mm=0; mm<length(p_Msg_prod0); mm++) {
        // Map sign-based index to magnitude-based index, restoring symmetry
        int mm_out = signed_to_unsigned_idx(mm, res_inputs);
        p_Msg_prod_combined(mm_out) += p_Msg_prod0(mm);
    }
    
    return p_Msg_prod_combined;
}

template <class Num_T> Vec<Num_T> lut_ldpc::unique(const Vec<Num_T>& x){
    if(length(x)<=1) return x;
    
    Vec<Num_T> y = x;
    sort(y);
    
    int ii=0;
    int len = length(y);
    while(ii < len-1){
        if( y(ii) == y(ii+1) )  y.del(ii);
        else ii++;
    }
    return y;
}
template  Vec<int> lut_ldpc::unique<int>(const Vec<int>& x);

inline double lut_ldpc::sig2snr(double rate, double sig){
    return -10*std::log10(2*rate*itpp::sqr(sig));
}

inline double lut_ldpc::snr2sig(double rate, double snr){
    return  itpp::pow10(-snr/20)/std::sqrt(2*rate);
}

vec lut_ldpc::sig2snr(double rate, const vec& sig){
    int len = sig.length();
    vec snr(len);
    for(int ii=0; ii<len; ii++)
        snr(ii) = sig2snr(rate, sig.get(ii));
    return snr;
}

vec lut_ldpc::snr2sig(double rate, const vec&  snr){
    int len = snr.length();
    vec sig(len);
    for(int ii=0; ii<len; ii++)
        sig(ii) = snr2sig(rate, snr.get(ii));
    return sig;
}

int lut_ldpc::quant_lin(double x, double delta, int N){
    it_assert_debug(mod(N,2) == 0, "quant_lin(): Number of quantization intervals must be even");
    int y = ceil_i(x/delta)+N/2-1;
    if(y<0) return 0;
    else if(y>N-1) return N-1;
    else return y;
}

int lut_ldpc::quant_nonlin(double x, const vec& boundaries){
    int idx = 0;
    for(int ii=0; ii<length(boundaries); ii++){
        if(x > boundaries(ii))
            idx++;
        else
            break;
    }
    return idx;
}

ivec lut_ldpc::quant_nonlin(const vec& x, const vec& boundaries){
    int len = length(x);
    ivec y(len);
    for(int ii=0; ii< len; ii++){
        y(ii) = quant_nonlin(x(ii), boundaries);
    }
    return y;
}

vec lut_ldpc::get_gaussian_pmf(double mu, double sig, int N, double delta){
    vec pmf(N);
    pmf(0) = 1-Qfunc(((-N/2.0+1)*delta-mu)/sig);
    for(int nn=1; nn<N-1; nn++){
        pmf(nn) = Qfunc( ((nn-N/2.0)*delta-mu)/sig) -
        Qfunc( ((nn+1-N/2.0)*delta-mu)/sig);
    }
    pmf(N-1) = Qfunc(((N/2.0-1)*delta-mu)/sig);
    return pmf/sum(pmf);
}


double lut_ldpc::rate_to_shannon_thr(double R){
    return 1.0/std::sqrt(pow2(2*R)-1);
}

double lut_ldpc::shannon_thr_to_rate(double sig){
    return .5*std::log2(1+ 1.0/sqr(sig));
}



inline double lut_ldpc::x_log2_y(double x, double y){
    if(x==0)   return 0;
    else if(x>0 && y>0) return x*std::log2(y);
    it_error("x_log2_y():Input invalid");
    return -1;
}


template <class Num_T> Vec<Num_T> lut_ldpc::fliplr(const Vec<Num_T>& x){
    int len = x.length();
    Vec<Num_T> y(len);
    for(int ii=0; ii<len; ii++ )    y(ii) = x(len-1-ii);
    return y;
}

template <class Num_T> Vec<Num_T> lut_ldpc::kron(const Vec<Num_T>& x, const Vec<Num_T>& y){
    int len_x = length(x);
    int len_y = length(y);
    Vec<Num_T> z(len_x*len_y);
    for(int xx=0; xx<len_x; xx++){
        for(int yy=0; yy<len_y; yy++){
            z(xx*len_y + yy) = x(xx)*y(yy);
        }
    }
    return z;
    
}

int lut_ldpc::signed_to_unsigned_idx(int idx, const ivec& inres){
    it_assert(idx < prod(inres) && idx >= 0, "lut_ldpc::signed_to_unsigned_idx(): Index out of range");
    
    int num_in = length(inres);
    ivec idx_in = zeros_i(num_in);
    
    int out_max = 2 * prod(inres/2);
    
    // Get input labels
    int idx_tmp = idx;
    for(int jj=0; jj<num_in; jj++){
        idx_in(jj) = idx_tmp % inres(jj);
        idx_tmp /= inres(jj);
    }
    
    // Get output labels and parity
    int parity = 0;
    int idx_out=0;
    int base = 1;
    for(int jj=0; jj<num_in; jj++){
        if(idx_in(jj) < inres(jj)/2 ){
            parity ^= 1;
            idx_out += base * ( inres(jj)/2 - 1 -  idx_in(jj) );
        }
        else{
            idx_out += base * ( idx_in(jj) - inres(jj)/2 );
        }
        base *= inres(jj)/2;
    }
    
    if(parity == 0)
        return out_max-1-idx_out;
    else
        return idx_out;
    
}

double lut_ldpc::quant_mi_sym(vec& p_out, ivec& Q_out, const vec& p_in, int Nq, bool sorted){
    int K = Nq;
    int M_in = length(p_in);
    
    // asserts
    it_assert(M_in%2 == 0, "quant_mi_sym(): Length of input pmf must be even");
    it_assert(K%2 == 0, "quant_mi_sym(): Number of output labels must be even");
    
    vec p_sorted;
    int M;
    ivec idx_in, idx_sorted;
    if(sorted == false){
        // ensure input is sorted strictly and eliminate duplicated entries
        // this symmetric approach takes care of occurances of LLR 0, equally splitting it among both
        // halves of the pmf domain
        p_sorted = sym_llr_sort_unique(p_in,idx_in,idx_sorted);
        M = p_sorted.length();
    }
    else{
        idx_in = linspace_fixed_step(0, M_in-1);
        idx_sorted = linspace_fixed_step(0, M_in-1);
        p_sorted = p_in;
        M = M_in;
    }
    
    
    // trivial case
    if(K>=M){
        Q_out.set_size(M_in);
        int outlabel = 0;
        for(int mm=0; mm<M_in/2; mm++){
            if(idx_sorted(mm) > outlabel)
                outlabel++;
            Q_out(idx_in(M_in-1-mm))   = K- 1 - outlabel;
            Q_out(idx_in(mm)) =  outlabel;
        }
        // Compute output pmf
        p_out = zeros(K);
        for(int mm=0; mm<M_in; mm++){
            p_out(Q_out(mm)) += p_in(mm);
        }
        return lut_ldpc::get_mi_bcpmf_sym(p_in);
    }
    
    // Precompute partial mutual information. This matrix is going to be upper diagonal because a >= ap
    mat g = zeros(M/2,M/2);
    for (int ap=0; ap<M/2; ap++) {
        double p_plus=0;
        double p_minus=0;
        for (int a=ap; a<M/2; a++) {
            p_plus  += p_sorted(M/2+a);
            p_minus += p_sorted(M/2-1-a);
            g(ap,a) = x_log2_y(p_plus, 2*p_plus/(p_plus+p_minus)) +  x_log2_y(p_minus, 2*p_minus/(p_plus+p_minus));
        }
    }
    
    
    // Compute state values S and transitions h
    mat S  = zeros(M/2,K/2);
    imat h = zeros_i(M/2,K/2);
    for(int a=0; a<(M-K)/2+1; a++){
        S(a,0) = g(0,a);
    }
    for(int zz=1; zz<K/2; zz++){
        for(int a=zz; a<zz+(M-K)/2+1; a++){
            S(a,zz) = -std::numeric_limits<double>::max();
            for(int ap=zz; ap<=a; ap++){
                double t = S(ap-1,zz-1) + g(ap,a);
                if (t > S(a,zz)) {
                    S(a,zz) = t;
                    h(a,zz) = ap;
                }
            }
        }
    }
    
    // compute optimal interval boundaries
    ivec astar = zeros_i(K/2+1);
    astar(K/2)=M/2;
    for(int kk=K/2-1; kk>0; kk--){
        astar(kk) = h(astar(kk+1)-1, kk);
    }
    
    // Build quantizer
    int outlabel = 0;
    Q_out.set_size(M_in);
    for(int mm=0; mm<M_in/2; mm++){
        if(idx_sorted(mm+M_in/2)-M/2 >= astar(outlabel+1))    outlabel++;
        Q_out(idx_in(M_in/2+mm))   = K/2     + outlabel;
        Q_out(idx_in(M_in/2-1-mm)) = K/2-1   -  outlabel;
    }
    // Compute output pmf
    p_out = zeros(K);
    for(int mm=0; mm<M_in; mm++){
        p_out(Q_out(mm)) += p_in(mm);
    }
    
    
    // Mutual information
    return S(M/2-1,K/2-1);
    
}

vec lut_ldpc::sym_llr_sort_unique(const vec& p_in, ivec& idx_in, ivec& idx_sorted, double llr_delta){
    int M_in = p_in.length();
    vec llr = log(p_in) - fliplr(log(p_in));
    idx_in = sort_index(llr);
    int mm=0;
    while(mm<M_in){
        int ll=1;
        while( mm+ll < M_in &&  llr(idx_in(mm)) == llr(idx_in(mm+ll)) ) ll++;
        idx_in.set_subvector(mm, idx_in.get(mm+sort_index(idx_in.get(mm, mm+ll-1)) ) );
        mm +=ll;
    }
    it_assert_debug(sum( to_ivec(idx_in+fliplr(idx_in) != M_in-1) )==0, "lut_ldpc::sym_llr_sort_unique(): Couldn't find symmetric permutation");
    
    idx_sorted.set_size(M_in/2);
    idx_sorted(0) = 0;
    double dupl = llr(idx_in(0));
    int dupl_idx= 0;
    int num_dupl = 0;
    for(int mm=1; mm<M_in/2; mm++){
        if(std::abs(llr(idx_in(mm)) - dupl) <= llr_delta){
            num_dupl++;
        }
        else{
            dupl_idx++;
        }
        idx_sorted(mm) = dupl_idx;
        dupl = llr(idx_in(mm));
    }
    
    idx_sorted = concat(idx_sorted, 2*max(idx_sorted)+1 - fliplr(idx_sorted));
    int M = M_in-2*num_dupl;
    vec p_sorted = zeros(M);
    for(int mm=0; mm<M_in; mm++){
        p_sorted(idx_sorted(mm)) += p_in(idx_in(mm));
    }
    return p_sorted;
}

double lut_ldpc::get_mi_bcpmf_sym(const vec& p){
    int K = length(p);
    it_assert(K>0 && mod(K, 2)==0, "get_mi_bcpmf_sym(): Input must have size >0 and an even number of elements");
    double mi = 0;
    for (int ii=0; ii<K/2; ii++) {
        mi +=     p(ii)*     std::log2(2*p(ii)/    ( p(ii)+p(K-1-ii)) ) +
        p(K-1-ii)* std::log2(2*p(K-1-ii)/( p(K-1-ii)+p(ii)) );
    }
    return mi;
}
