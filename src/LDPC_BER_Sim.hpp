/*!
 * \file
 * \brief
 * \author Michael Meidlinger
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2016 Michael Meidlinger - All Rights Reserved
 *
 */


#ifndef LDPC_BER_Sim_hpp
#define LDPC_BER_Sim_hpp

#include "LDPC_DE.hpp"
#include "LDPC_Code_LUT.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <ctime>
#include <stdio.h>
#include <sys/stat.h>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/optional.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iomanip>

//! Git version of program
extern const char *gitversion;

namespace itpp{

    
class LDPC_BER_Sim_Results{
public:
    LDPC_BER_Sim_Results(){}
    LDPC_BER_Sim_Results(int nvar, int nchk): ldpc_nvar(nvar), ldpc_nchk(nchk)
        {ldpc_code_rate = 1.0- (double)nchk/nvar;}
    //! Write Params to .m file
    void save();
    void write_itfile(const std::string& filename);
    void add_snr_point(double snr, const BERC& berc_data, const BERC& berc_uncoded, const BLERC& ferc);
    void save_runtime(double runtime){this->runtime = runtime;}
protected:
    
    int     ldpc_nvar;
    int     ldpc_nchk;
    double  ldpc_code_rate;
    double  runtime;                //!< Runtime of the simulation in seconds

    
    vec             sim_SNRdB;
    Vec<int64_t>    sim_Nframes;
    Vec<int64_t>    sim_Ndatabits;
    Vec<int64_t>    sim_frame_errors;
    Vec<int64_t>    sim_data_bit_errors;
    Vec<int64_t>    sim_uncoded_bit_errors;
    
    
};
   
/*!
 \brief Base class for BER Simulations Parameters
*/
class LDPC_BER_Sim{
    
public:
    //! Default Constructor
    LDPC_BER_Sim(): decoder_set(false), encoder_set(false), sim_finished(false){}
    
    //! Constructor from param file name. The param file name is
    LDPC_BER_Sim(const boost::filesystem::path& params_file_path, const boost::filesystem::path& base_dir_path);
    
    //! Virtual Destructor
    virtual ~LDPC_BER_Sim() {};

    //! Run Simulation
    void run();
    //! Load LDPC Code
    virtual void load();
    //! Save simulation files
    void save();
    
    //! Simulate a specific SNR point.
    bool sim_snr_point( double snr);
    
    //! Generate filename
    virtual std::string gen_filename() const;

    //! Extend the custom file name appended at the end of the results filename
    void append_custom_name(std::string& custom_name_ext);
    
    
    // Simulation Parameters
    int64_t         Nframes;                //!< Number of frames to simulate
    int             Nfers;                  //!< Number of frame errors after which the simulation will be terminated
    vec             SNRdB;                  //!< SNR Points to simulate
    double          fer_min;                //!< Break once the frame error rate is below this threshold
    double          ber_min;                //!< Once a bit error rate of ber_min is reached, points with higher SNR are not processed any more
    bool            zero_codeword;          //!< If true, rather than a randomly generated encoded bitstring, the all zero word will be transmitted
    int             rand_seed;              //!< Random seed of simulation. Usually read in from argv[]
    int             rand_seed_offset;       //!< Offset added to random seed of simulation
    int             save_codec;             //!< For which value of the random seed the codec is saved. This is to avoid race conditions in case of multiple instances of the program wanting to write the codec at the same time. Set to a negative number to not save at all
    int             max_iter;               //!< Maximum Number of LDPC Iterations
    bool            parity_check_iter;      //!< Wether parity checks are performed after each decoding iteration
    std::string     parity_filename;        //!< Filename of Parity Check Matrix without .alist extension (alist format assumed)
    std::string     custom_name;            //!< This string will be attached to the auto generated filename for the resuls file
    std::string     results_prefix;         //!< Any results filename and folder will start with this prefix
    bool            save_permuted;          //!< If the simulation constructs a Generator it usually has to permute columns of the parity check Matrix. If this flag is true, the permuted versions of both Generator and Parity are stored in codes_dir, potentially overwriting the originals Parity .alist file.
    
    // Search paths for the program
    std::string     results_dir;
    std::string     codes_dir;
  
    boost::filesystem::path results_path;
    boost::filesystem::path codes_path;
    boost::filesystem::path params_file_path;
    std::string codec_filename;
    
    

    
    
    // Parameters of the code
    int codeword_length;
    int dataword_length;
    
    bool decoder_set;
    bool encoder_set;
    bool sim_finished;
    
protected:
    // LDPC code implementation.
    Channel_Code* C;
    LDPC_BER_Sim_Results results;
    

    
private:
    /*!
        \brief Internal LLR resolution
    C.f. http://itpp.sourceforge.net/4.3.1/classitpp_1_1LLR__calc__unit.html#a531b0a4eca593e439a39ba9de42ee68c for details.
     */
    int llr_calc_d1;
    //! Resolution of JacLog LUT. Set to 0 to emulate a Min-Sum Decoder
    int llr_calc_d2;
    int llr_calc_d3;
    int llr_calc_d4;
};
    
    
class LDPC_BER_Sim_LUT: public LDPC_BER_Sim{
public:
    LDPC_BER_Sim_LUT(const boost::filesystem::path& params_file_path, const boost::filesystem::path& base_dir_path);
    
    
    int tree_autogen_mode;                  //!< Wether the trees are generated automatically or have to be specified manually
    bvec reuse_lut;                          //!< This vector is 1 at position ii if at iteration ii, the LUT of iteration ii-1 is reused, ii=1,...
    int Nq_Cha;                             //!< Resolution of channel input LLRs
    ivec Nq_Msg;                            //!< Quantizer ii has input resolution Nq_Msg(ii) and output resolution Nq_Msg(ii+1), ii=1,...
    bool min_lut;                            //!< If true, the check node update is based on a minimum operation on sorted labels rather than a LUT
    boost::optional<double> design_thr;      //!< For which noise standard deviation (Assuming BIAWGN  hannel) the LUT's should be designed.
    boost::optional<double> design_SNRdB;    //!< For which (Eb/N0)_dB (Assuming BIAWGN  hannel) the LUT's should be designed.
    
    /*!
     If \c tree_auto_gen == false, then the trees have to be specified manually here, one for each possible node degree
    */

    //! Generate filename
    virtual std::string gen_filename() const;
    
    std::string tree_mode;
    std::string trees_filename;
    std::string     trees_dir;
    boost::filesystem::path trees_path;
    
    

    //! Check if parameters match
    void check_params();
    virtual void load();
    /*! \brief Design Quantizers for LUT decoders by means of density evolution
     */
    LDPC_Code_LUT* design_lut_codec(const LDPC_Parity *H, LDPC_Generator *G) const;
};

}
#endif
