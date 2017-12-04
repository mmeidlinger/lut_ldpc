/*!
 * \file
 * \brief Program to convert .ens files to .deg files
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
