/*!
 * \file
 * \brief Program to convert .ens files to .deg files
 * \author Michael Meidlinger
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2016 Michael Meidlinger - All Rights Reserved
 *
 */

#include <itpp/itbase.h>
#include "LDPC_Ensemble.hpp"

using namespace lut_ldpc;
using namespace itpp;
using namespace std;

int main(int argc, char **argv){
    // Parse file name
    if(argc!=3){
        cout << "Usage: ens2deg infile.ens outfile.deg" << endl;
        return EXIT_FAILURE;
    }
    
    string ensfilename(argv[1]);
    string degfilename(argv[2]);
    // Load ensemble
    LDPC_Ensemble ens(ensfilename);
    ens.export_deg(degfilename);
    return EXIT_SUCCESS;
}
