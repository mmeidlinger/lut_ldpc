/*!
 * \file
 * \brief Program to generate an ensemble from a parity check matrix.
 * \author Michael Meidlinger
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2016 Michael Meidlinger - All Rights Reserved
 *
 */

#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include "LDPC_DE.hpp"
#include "LDPC_Code_LUT.hpp"
#include <boost/format.hpp>

using namespace itpp;
using namespace std;
int main(int argc, char **argv){

    // Parse input name
    if(argc!=3){
        cout << "Usage: gen_ensemble infile.alist outfile.ens" << endl;
        return EXIT_FAILURE;
    }
    
    string parity_filename(argv[1]);
    string ens_filename(argv[2]);
    

    LDPC_Parity H(parity_filename);

    cout << "Calculating LDPC degree distributions ..." << endl;
    LDPC_Ensemble ens = get_empirical_ensemble(H);
    cout << "Done. Rate of Ensemble = " << ens.get_rate() << endl;
    cout << ens;
    ens.write(ens_filename);

    return EXIT_SUCCESS;
}
