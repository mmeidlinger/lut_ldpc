//
//  LDPC_Ensemble.cpp
//  LDPC_LUT
//
//  Created by Michael Meidlinger on 28.11.17.
//

#include "LDPC_Ensemble.hpp"
#include "TextTable.hpp"
#include "common.hpp"


using namespace lut_ldpc;
using namespace itpp;
using namespace std;

// ----------------------------------------------------------------------
// LDPC_Ensemble
// ----------------------------------------------------------------------

const double LDPC_Ensemble::pmass_tolerance = 1e-2;

LDPC_Ensemble::LDPC_Ensemble(): init_flag(false){}

LDPC_Ensemble::LDPC_Ensemble(const vec& l, const vec& r): init_flag(false){
    set_chk_degree_dist(r);
    set_var_degree_dist(l);
    check_consistency();
    init_flag = true;
}

void LDPC_Ensemble::set_chk_degree_dist(const vec& r)
{
    dc_act = 0;
    for(int ii=0; ii< length(r); ii++)
        if(r(ii) > 0)   dc_act++;
    
    rho = zeros(dc_act);
    degree_rho = ivec(dc_act);
    
    int idx=0;
    for(int ii=0; ii< length(r); ii++){
        if(r(ii) > 0){
            degree_rho(idx) = ii+1;
            rho(idx) = r(ii);
            idx++;
        }
    }
    if(init_flag) check_consistency();
}

void LDPC_Ensemble::set_var_degree_dist(const vec& l)
{
    dv_act = 0;
    for(int ii=0; ii< length(l); ii++)
        if(l(ii) > 0)   dv_act++;
    
    lam = zeros(dv_act);
    degree_lam = ivec(dv_act);
    
    int idx=0;
    for(int ii=0; ii< length(l); ii++){
        if(l(ii) > 0){
            degree_lam(idx) = ii+1;
            lam(idx) = l(ii);
            idx++;
        }
    }
    if(init_flag) check_consistency();
}

inline void LDPC_Ensemble::check_consistency(){
    
    // Check if degree distributions are nonnegative
    bool lam_is_pos, rho_is_pos;
    lam_is_pos = !(bool)sum(to_ivec(lam<0));
    rho_is_pos = !(bool)sum(to_ivec(rho<0));
    it_assert( lam_is_pos && rho_is_pos, "LDPC_Ensemble::check_consistency(): Degree distributions need to be positive!" );
    
    // Check if degrees are uniqe
    it_assert( degree_lam.length() == unique(degree_lam).length() &&
              degree_rho.length() == unique(degree_rho).length() ,
              "LDPC_Ensemble::check_consistency(): Degree vectors need to have unique entries" );
    
    // Check if dimensions match
    it_assert( dv_act == length(lam) &&
              dv_act == length(degree_lam) &&
              dc_act == length(rho) &&
              dc_act == length(degree_rho),
              "LDPC_Ensemble::check_consistency(): Inconsistent dimension");
    
    // Check if degrees  are non_negative
    lam_is_pos = !(bool)sum(to_ivec(degree_lam<0));
    rho_is_pos = !(bool)sum(to_ivec(degree_rho<0));
    it_assert( lam_is_pos && rho_is_pos, "LDPC_Ensemble::check_consistency(): Degree distributions need to be positive!" );
    
    
    // Check if degree distributions sum to one
    double sum_lam = sum(lam);
    double sum_rho = sum(rho);
    it_assert( std::abs(1.0-sum_lam) < pmass_tolerance ||
              std::abs(1.0-sum_rho) < pmass_tolerance ,
              "LDPC_Ensemble::check_consistency(): Degree distributions do not sum to one within tolerance!" );
    
    // Normalize
    lam = lam/sum_lam;
    rho=rho/sum_rho;
    
    // Check if Coderate is positive
    it_assert(get_rate()>0, "LDPC_Ensemble::check_consistency(): Coderate is neagtive!");
}

LDPC_Ensemble::LDPC_Ensemble(const ivec& dl, const vec& l, const ivec& dr, const vec& r): init_flag(false){
    dv_act = dl.length();
    dc_act = dr.length();
    
    it_assert( sum(to_ivec(dl < 1)) == 0 && sum(to_ivec(dr < 1)) == 0, "LDPC_Ensemble::LDPC_Ensemble(): Degrees must be larger than 0" );
    it_assert( dv_act == length(l) && dc_act == length(r), "LDPC_Ensemble::LDPC_Ensemble(): Input dimension mismatch");
    
    degree_lam = dl;
    degree_rho = dr;
    rho = r;
    lam = l;
    
    check_consistency();
    init_flag = true;
}

LDPC_Ensemble::LDPC_Ensemble(const std::string& filename): init_flag(false)
{
    read(filename);
}

void LDPC_Ensemble::read(const std::string& filename)
{
    std::ifstream file;
    std::string line;
    std::stringstream ss;
        
    file.open(filename.c_str());
    it_assert(file.is_open(),
              "LDPC_Ensemble::read(): Could not open file \""
              << filename << "\" for reading");
    
    // parse number of active degrees
    getline(file, line);
    ss << line;
    ss >> dv_act >> dc_act;
    it_assert(!ss.fail(),
              "LDPC_Ensemble::read(): Wrong active degree data!");
    it_assert( (dv_act > 0) && (dc_act > 0),
              "LDPC_Ensemble::read(): Wrong active degree data!");
    
    ss.seekg(0, std::ios::end);
    ss.clear();
    
    
    // parse variable node distribution
    degree_lam.set_size(dv_act);
    degree_lam.clear();
    lam.set_size(dv_act);
    lam.clear();
    getline(file, line); //variable node degrees
    ss << line;
    for (int ii = 0; ii < dv_act; ii++) {
        ss >> degree_lam(ii);
        it_assert(!ss.fail(),
                  "LDPC_Ensemble::read(): Wrong active degree data! (degree_lam("
                  << ii << "))");
        it_assert((degree_lam(ii) >= 1),
                  "LDPC_Ensemble::read(): Wrong active degree data! (degree_lam("
                  << ii << "))");
    }
    ss.seekg(0, std::ios::end);
    ss.clear();
    getline(file, line); //variable node distribution
    ss << line;
    for (int ii = 0; ii < dv_act; ii++) {
        ss >> lam(ii);
        it_assert(!ss.fail(),
                  "LDPC_Ensemble::read(): Wrong active degree data! (lam("
                  << ii << "))");
        it_assert((lam(ii) > 0) && (lam(ii) <= 1),
                  "LDPC_Ensemble::read(): Wrong active degree data! (lam("
                  << ii << "))");
    }
    ss.seekg(0, std::ios::end);
    ss.clear();
    
    // parse check node distribution
    degree_rho.set_size(dc_act);
    degree_rho.clear();
    rho.set_size(dc_act);
    rho.clear();
    getline(file, line); //check node degrees
    ss << line;
    for (int ii = 0; ii < dc_act; ii++) {
        ss >> degree_rho(ii);
        it_assert(!ss.fail(),
                  "LDPC_Ensemble::read(): Wrong active degree data! (degree_rho("
                  << ii << "))");
        it_assert((degree_rho(ii) >= 1),
                  "LDPC_Ensemble::read(): Wrong active degree data! (degree_rho("
                  << ii << "))");
    }
    ss.seekg(0, std::ios::end);
    ss.clear();
    getline(file, line); //check node distribution
    ss << line;
    for (int ii = 0; ii < dc_act; ii++) {
        ss >> rho(ii);
        it_assert(!ss.fail(),
                  "LDPC_Ensemble::read(): Wrong active degree data! (rho("
                  << ii << "))");
        it_assert((rho(ii) > 0) && (rho(ii) <= 1),
                  "LDPC_Ensemble::read(): Wrong active degree data! (rho("
                  << ii << "))");
    }
    ss.seekg(0, std::ios::end);
    ss.clear();
    
    file.close();
    
    check_consistency();
    init_flag = true;
    return;
}

void LDPC_Ensemble::write(const std::string& filename) const
{
    std::ofstream file;
    
    
    file.open(filename.c_str());
    it_assert(file.is_open(),
              "LDPC_Ensemble::write(): Could not open file \""
              << filename << "\" for writing");
    
    // write number of active degrees
    file <<   dv_act << " " << dc_act << endl;
    
    // write variable node degrees and distribution
    for (int ii = 0; ii < dv_act-1; ii++) {
        file << degree_lam(ii) << " ";
    }
    file << degree_lam(dv_act-1)<< endl;
    
    for (int ii = 0; ii < dv_act-1; ii++) {
        file << lam(ii) << " ";
    }
    file << lam(dv_act-1)<< endl;
    
    // write variable node degrees and distribution
    for (int ii = 0; ii < dc_act-1; ii++) {
        file << degree_rho(ii) << " ";
    }
    file << degree_rho(dc_act-1)<< endl;
    
    for (int ii = 0; ii < dc_act-1; ii++) {
        file << rho(ii) << " ";
    }
    file << rho(dc_act-1)<< endl;
    
    file.close();
    return;
}

void LDPC_Ensemble::export_deg(const std::string& filename) const{
    std::ofstream file;
    
    
    file.open(filename.c_str());
    it_assert(file.is_open(),
              "LDPC_Ensemble::export_deg(): Could not open file \""
              << filename << "\" for writing");
    
    // write number of active variable node degrees
    file << dv_act << endl;
    
    
    // write variable node degrees
    for (int ii = 0; ii < dv_act-1; ii++) {
        file << degree_lam(ii) << " ";
    }
    file << degree_lam(dv_act-1)<< endl;
    
    // write variable node degree distribution (node perspective)
    vec Lam = sget_Lam();
    it_assert(Lam.size()==dv_act, "LDPC_Ensemble::export_deg(): Inconsistent ensemble");
    for (int ii = 0; ii < dv_act-1; ii++) {
        file << Lam(ii) << " ";
    }
    file << Lam(dv_act-1)<< endl;
    
    file.close();
    return;
}

double LDPC_Ensemble::get_rate() const {
    return 1 -  sum(elem_div(rho, to_vec(degree_rho))) / sum(elem_div(lam, to_vec(degree_lam)));
}

vec LDPC_Ensemble::sget_lam()const  {return lam;}

int LDPC_Ensemble::sget_lam(vec& l, ivec& dl)const{
    l = lam;
    dl = degree_lam;
    return dv_act;
}

vec LDPC_Ensemble::sget_Lam()const {
    vec Lam(dv_act);
    Lam = elem_div(lam, to_vec(degree_lam));
    return Lam/sum(Lam);
}

vec LDPC_Ensemble::sget_rho()const{ return rho; }

int LDPC_Ensemble::sget_rho(vec& r, ivec& dr)const{
    r = rho;
    dr = degree_rho;
    return dc_act;
}

vec LDPC_Ensemble::sget_Rho()const {
    vec Rho(dc_act);
    Rho = elem_div(rho, to_vec(degree_rho));
    return Rho/sum(Rho);
}

ivec LDPC_Ensemble::sget_degree_rho() const {return degree_rho;}
ivec LDPC_Ensemble::sget_degree_lam() const {return degree_lam;}

int LDPC_Ensemble::get_dv_act() const   {return dv_act;}
int LDPC_Ensemble::get_dc_act() const   {return dc_act;}


vec LDPC_Ensemble::get_chk_degree_dist() const {
    vec r = zeros(max(degree_rho));
    for(int ii=0; ii<dc_act; ii++) r(degree_rho(ii)-1) = rho(ii);
    return r;
}

vec LDPC_Ensemble::get_var_degree_dist() const {
    vec l = zeros(max(degree_lam));
    for(int ii=0; ii<dv_act; ii++) l(degree_lam(ii)-1) = lam(ii);
    return l;
}

void LDPC_Ensemble::sset_rho(vec r){
    it_assert_debug(r.length() == dc_act, "LDPC_Ensemble::sset_rho(): Active number of degrees does not match");
    rho = r;
    if(init_flag) check_consistency();
}

void LDPC_Ensemble::sset_lam(vec l){
    it_assert_debug(l.length() == dv_act, "LDPC_Ensemble::sset_rho(): Active number of degrees does not match");
    lam = l;
    if(init_flag) check_consistency();
}

double LDPC_Ensemble::get_lam_of_degree(int d) const {
    
    for(int ii=0; ii<dv_act; ii++){
        if( d == degree_lam(ii) )   return lam(ii);
    }
    return 0;
}

LDPC_Ensemble lut_ldpc::get_empirical_ensemble(const LDPC_Parity& H){
    static const int max_degree = 200;
    
    
    
    int nvar = H.get_nvar();
    int nchk = H.get_ncheck();
    
    // Iterate over all elements. This can be done smarter and faster, however, for that we would need to access internals of the LDPC_Parity class
    
    
    
    ivec col_sum = H.get_colsum();
    ivec row_sum = H.get_rowsum();
    
    vec var_edge_deg = zeros(max_degree);
    vec chk_edge_deg = zeros(max_degree);
    
    for(int nn=0; nn< nvar; nn++){
        it_assert_debug(col_sum(nn) <= max_degree, "get_empirical_ensemble(): Maximum degree exceeded");
        it_assert_debug(col_sum(nn) > 0, "get_empirical_ensemble(): Minimum degree is 1");
        var_edge_deg(col_sum(nn)-1)+= col_sum(nn);
    }
    for(int mm=0; mm< nchk; mm++){
        it_assert_debug(row_sum(mm) <= max_degree, "get_empirical_ensemble(): Maximum degree exceeded");
        it_assert_debug(row_sum(mm) > 0, "get_empirical_ensemble(): Minimum degree is 1");
        chk_edge_deg(row_sum(mm)-1)+= row_sum(mm);
    }
    
    var_edge_deg = var_edge_deg/sum(var_edge_deg);
    chk_edge_deg = chk_edge_deg/sum(chk_edge_deg);
    return LDPC_Ensemble(var_edge_deg, chk_edge_deg);
}
