/*!
 * \file
 * \brief Program to generate an ensemble from a parity check matri in alist format
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
#include <itpp/itcomm.h>
#include "LDPC_Ensemble.hpp"
#include <boost/format.hpp>

using namespace itpp;
using namespace std;
using namespace lut_ldpc;

int main(int argc, char **argv){

    // Parse input name
    if(argc!=3){
        cout << "Usage: alist2ensemble infile.alist outfile.ens" << endl;
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
