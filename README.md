# Introduction

LUT LDPC is a collection of software tools to design and test LDPC decoders based on discrete message passing decoding
using lookup tables (LUTs), referred to as LUT decoders, cf. [[1]](#literature). It is mainly written in C++ and relies  on the [IT++](http://itpp.sourceforge.net/) for abstracting basic linear algebra and signal processing operations.
Consequently, the LUT decoders can easily be integrated into more complex communication systems including concatenated coding and/or modulation.

# Installation

## Requirements
The [boost C++ libraries](http://www.boost.org/) are used for file system abstraction and program I/O.
[IT++](http://itpp.sourceforge.net/) is included via a submodlue because we had to do some patching to get access to internals of its classes.
In order to successfully compile with IT++, you require its dependencies that are FFTW3, LAPACK an BLAS, cf. [http://itpp.sourceforge.net/4.3.1/](http://itpp.sourceforge.net/4.3.1/). Often, those packages are included with your OS already, but usually as dynammically linked libraries.

The program has been tested on MacOSX and Linux, but with minor modification, will most likely also run on Windows.
The following instructions  refer to a Linux install.

## Installing Dependencies
On modern Linux distributions you can fetch all neccesary software to build LUT LDPC and required libraries using a package manager. E.g., on Ubuntu you can
install them via
```bash
$ sudo apt-get install git build-essential cmake libboost-all-dev libfftw3-dev liblapack-dev libblas-dev
```

## Cloning and Building
We favoured static over dynamic linking for portability, as we ran the compiled binaries on an inhomogeneous cluster of Linux hosts.
Dynamically linked versions could be obtained by changing the [Makefile](./Makefile) appropriately.

To download and install the software open a terminal and enter
```bash
git clone --recursive https://github.com/mmeidlinger/lut_ldpc.git
cd lut_ldpc
make
make install
```
This will compile and link the program as well as the patched IT++ library and install the programs locally within the `lut_ldpc` directory into the `bin` directory.
By default, the libraries that are usually not included within your OS (boost and itpp) are linked statically to the binary for reasons of portability.
So `make` is shorthand for  `make BUILD_TYPE=Release LINK_TYPE=static INSTALLDIR=bin`.


# Usage

LUT-LDPC consists of several programs and scripts, whose usage we will discuss in what follows. In general, we want to emphasize that LUT-LDPC is a research driven toolset. As such, its functionality should be used more like a library than via the included programs. Hence,
the ways to interact with the programs are limited and applying the software to specific and complex scenarios requires the user to write and compile their own programs. To assist users with that, LUT-LDPC features [Doxygen](http://www.doxygen.org/) inline source code documentation.

## Designing LUT Decoders and testing Bit Error rate performance
### Running the simulation
`ber_sim` is capable of designing decoders and conducting bit error rate (BER) Monte Carlo simulations.
It takes the following parameters, which can be displayed by
```bash
$ bin/ber_sim -h
OPTIONS:
    -b [ --basedir ] arg (=/absolute/path/to/lut_ldpc)  paths in params files are relative to this directory. Default: current direcroy
    -c [ --custom-name ] arg                            append this string at the end of the results file name
    -h [ --help ]                                       produce help message
    -p [ --params ] arg                                 input parameter file
    -s [ --seed ] arg (=0)                              random seed
```
Most importantly, the `-p` option specifies the parameter file to configure the simulation. An example of such a files is given below, which is a stripped down version of
[`params/ber.ini.irregular.example`](params/ber.ini.irregular.example):
```
[Sim]
   SNRdB    = 0:.5:4
   Nframes  = 1e2

[LDPC]
   parity_filename = rate0.50_dv02-17_dc08-09_lut_q4_N500

[LUT]
   max_iter = 50
   design_thr = 0.88
   qbits_channel = 4
   qbits_message_uniform = 4
```
According to these settings, the `ber_sim` constructs a LUT decoder for an AWGN channel with noise standard deviation 0.88 (`design_thr = 0.88`). The channel LLRs are quantized with 4 bits resolution (`qbits_channel = 4`)  and the message resolution is 4 bits throughout the decoding process (`qbits_message_uniform = 4`).
For the LDPC code, the parity check matrix `codes/rate0.50_dv02-17_dc08-09_lut_q4_N500.alist` is used (`parity_filename = rate0.50_dv02-17_dc08-09_lut_q4_N500`)
and for each SNR point (`SNRdB    = 0:.5:4`) , at most 100 frames are simulated (`Nframes  = 1e2`). The results as well as the decoder object of the simulation are automatically saved into the `results` directory, where the name is derived from the simulation setting ( The `custom_name =` option can be used to apped a string to the auto generated results name).

So running
```bash
$ bin/ber_sim -p params/ber.ini.irregular.example
```
Produces
```
 results/
    RES_N500_R0.5_maxIter50_zcw0_frames100_minLUT/
        ber.ini.irregular.example
        lut_codec.it
        RES_N500_R0.5_maxIter50_zcw0_frames100_minLUT_rseed0000.it
```
As we can see, the folder `RES_N500_R0.5_maxIter50_zcw0_frames100_minLUT` is generated containing the original parameter file, the decoder object (`lut_codec.it`) as well as the results file `RES_N500_R0.5_maxIter50_zcw0_frames100_minLUT_rseed0000.it`.

### Evaluating the results
The results can now be visualized using the MATLAB script  [`scripts/analyze_results.m`](scripts/analyze_results.m)
```matlab
K>> addpath scripts;
K>> analyze_results('results, {RES_N500_R0.5_maxIter50_zcw0_frames100_minLUT_rseed0000.it}')
```

### Comparing to other decoders
Have a look at [`params/ber.ini.irregular.example`](params/ber.ini.irregular.example) for different configuration options. `ber_sim` also supports conventional belief propagation (BP) decoding and min-sum decoding. We patched IT++ to support decoding with finite bit width to compare low resolution BP and LUT decoding.

## Density evolution
### Running simulations
`de_sim` is the program to run density evolution simulation. Except for the `-h` option,  in only takes a parameter file  via the `-p`  option:
```bash
$ bin/ber_sim -p params/de_sim.ini.example
```
To simulate [`params/de.ini.example`](params/de.ini.example)
```
[Sim]
thr_min = 1e-7
thr_prec = 1e-5
Pe_max = 1e-10
maxiter_de = 2000
max_ni_de_iters = 1
LLR_max = 25.0
results_name = results/rate0.50_dv02-17_dc08-09_lut_q4_minLUT1_de-auto-bin-balanced_example.txt
ensemble_filename = ensembles/rate0.50_dv02-17_dc08-09_lut_q4.ens

[LUT]
min_lut = true
qbits = 4 4
tree_mode = auto_bin_balanced
irregular_design_strategy = joint_root
```
Produces the output `results/rate0.50_dv02-17_dc08-09_lut_q4_minLUT1_de-auto-bin-balanced_example.txt` containing
```
==== DE Threshold for ensemble file ensembles/rate0.50_dv02-17_dc08-09_lut_q4.ens (Rate = 0.5, BI-AWGN channel)
Active Variable node degrees: [2 3 9 17]
pmf of Variable node edges: [0.138045 0.401038 0.026586 0.434331]
Active Check node degrees: [8 9]
pmf of Check node edges: [0.323376 0.676624]
-- SIMULATION PARAMETERS  Search Window = [1e-07, 1]
Threshold precision = 1e-05
Convergence error probability = 1e-10
Maximum Number of message passing iterations = [2000]
MinLut Algorithm used = 1
LUT Tree design mode = auto_bin_balanced
LUT table design mode = joint_root
LUT reuse iter vec = [0]
Non improving iterations tolerated before terminating = 1
Resolutions [channel bits, message bits; ...] = [[4 4]]
Program git version = 2809f10085499f057cb4e5f38f24afc0afe7c406
Bisection iterations until convergence = [20]
Stable lam2 degrees at thresholds = [0.131418]
Threshold(s) found = [0.929193]
Eb/N0 corresponding to thresholds = [0.637884]

```


## Generating codes
In [[2]](#literature), we found out that for irregular codes under LUT decoding, degree distributions must be optimized. To generate codes from optimized degree distributions,
this repository contains a copy of the PEG [[4]](#literature) program which is freely available at http://www.inference.org.uk/mackay/PEG_ECC.html
The copy resides in the `peg` subdirectory and must be compiled separately
```bash
$ cd peg
$ make
```
Since the peg program expects degree distributions in a slightly different format and doesn't use the `.alist` format for outputting parity check matrices, we wrote the [peg.sh](scripts/peg.sh) script
to directly convert ensembles to `.alist` files. The script makes use of the programs [`ens2deg`](prog/ens2deg.cpp) and [`dat2alist`](prog/dat2alist.cpp) to convert intputs and outputs to the peg program. E.g., to create a parity check matrix for a length 1000 rate 1/2 code from the degree distribution `ensembles/rate0.50_dv02-17_dc08-09_lut_q4.ens`
and save it under `codes/rate0.50_dv02-17_dc08-09_lut_q4_N1000.alist`,
```bash
$ scripts/peg.sh 500 1000 codes/rate0.50_dv02-17_dc08-09_lut_q4_N1000.alist ensembles/rate0.50_dv02-17_dc08-09_lut_q4.ens
```

## Optimizing the reuse pattern of LUTs
In general, every iteration of a LUT decoder implements different LUTs. LUTs can be reused for more than one iteration, however, the pattern of reuse must be carefully designed.
We provide the [`reuse_vec_opt`](prog/reuse_vec_opt.cpp) program to do this:
```bash
$ bin/reuse_vec_opt --help
    Called program via
    reuse_vec_opt --help
    OPTIONS:
        -m [ --min-approx ]                         Approximate Check Node Updates
        --quant-bits-msg arg (=4)                   Number of quantization bits for messages
        --quant-bits-cha arg (=4)                   Number of quantization bits for channel outputs
        -t [ --threshold ] arg                      Noise value to run DE. If not provided, found by bisction bisec_search
        -e [ --ensemble ] arg                       Filename for initial ensemble
        -i [ --iterations ] arg (=100)              Number of Message passing iterations
        -d [ --degree-dist ] arg                    Degree ditribution is the form, "VN_degrees / VN_probabilities / CN_degrees / CN_probabilities "
        -s [ --scale-down ] arg (=0.995)            Scale down threshold by this value if an updated ensemble fails to converge
        -p [ --pmax ] arg (=1e-11)                  Convergence error probability
        -r [ --reuse-stages ] arg                   Number of distinct LUT stages
        -v [ --reuse-vec ] arg                      Provide an initial reuse vector
        --lut-table-design arg (=joint_root)        Strategy for LUT table creation
        --lut-tree-design arg (=auto_bin_balanced)  Strategy for LUT table creation
        -h [ --help ]                               produce help message
```
Warning: `reuse_vec_opt` creates as many threads as iterations specified via the `-i` options, i.e., it is quite demanding on your CPU if executed on a client machine.
E.g. to optimize the reuse pattern for a code with degree distribution as specified in `ensembles/rate0.50_dv02-17_dc08-09_lut_q4.ens`,
```bash
$ bin/reuse_vec_opt -e ensembles/rate0.50_dv02-17_dc08-09_lut_q4.ens -t0.89 -i100  -p1e-11 -s0.999 --reuse-stages 25
```
Here, the program initially allocates 100 unique LUT stages for every iteration and tries to reduce the number of unique LUTs down to 25.


# Codes, Ensembles and Trees
C.f.  [codes](codes/README.md),  [ensembles](ensembles/README.md) and  [trees](trees/README.md).


# Creating VHDL Code for an Unrolled Decoder
Using the MATLAB scripts of the [LUT-LDPC-VHDL](lut_ldpc_vhdl/README.md) submodule, the VHDL source code for an unrolled decoder (cf. [[3-4]](#literature))
can be generated based on the decoders exported by LUT LDPC C++ program.
In order to generate a decoder object and some input-output pairs run

```bash
$ bin/ber_sim -p params/ber.ini.regular.example > stimuli.txt
```
This will simulate 20 frame per SNR value and create the decoder object  `results/RES_N2048_R0.841309_maxIter8_zcw0_frames20_minLUT/lut_codec.it`. Furthermore,  pairs of decoder inputs and corresponding decoder outputs are written to the text file `stimuli.txt`.

The decoder object  `lut_codec.it` can then be used as an input to generate corresponding VHDL, while the pairs in `stimuli.txt` can be used to verify the correctnes of the resulting VHDL against the simulation model.

For more details, cf.  the [documentation of LUT-LDPC-VHDL](lut_ldpc_vhdl/README.md)

# Writing your own programs
The [Makefile](./Makefile) is configured to compile one executable per source file in the `prog` directory and link it to all object files of LUT LDPC. Try adding [this](trees/README.md) example as `prog/tree_example.cpp` and rebuild and install using `make && make install`. This should give you the  program `bin/tree_example`.



# Referencing
If you use this software for your academic research, please consider referencing our original contributions [[1,2,3,4]](#literature)
# Literature
[[1] M. Meidlinger, A. Balatsoukas-Stimming, A. Burg, and G. Matz, “Quantized message passing for LDPC codes,” in Proc. 49th Asilomar Conf. Signals, Systems and Computers, Pacific Grove, CA, USA, Nov. 2015.](http://ieeexplore.ieee.org/document/7421419/)

[[2] M. Meidlinger and G. Matz, “On irregular LDPC codes with quantized message passing decoding,” in Proc. IEEE SPAWC 2017, Sapporo, Japan, Jul. 2017.
](http://ieeexplore.ieee.org/document/8227780/)

[[3]  A. Balatsoukas-Stimming, M. Meidlinger, R. Ghanaatian, G. Matz, and A. Burg, “A fully-unrolled LDPC decoder based on quantized message passing,” in Proc. SiPS 2015, Hang Zhou, China, 10 2015.
](http://ieeexplore.ieee.org/abstract/document/7345024/)

[[4] R. Ghanaatian, A. Balatsoukas-Stimming, C. Mu ̈ller, M. Meidlinger, G. Matz, A. Teman, and A. Burg, “A 588 Gbps LDPC decoder based on finite-alphabet message passing,” IEEE Trans. VLSI Systems, vol. 26, no. 2, pp. 329–340, 2 2018.
](http://ieeexplore.ieee.org/document/8113527/)

[[5] X.-Y. Hu, E. Eleftheriou, and D. Arnold, “Regular and irregular progressive edge-growth tanner graphs,” IEEE Trans. Information Theory, vol. 51, no. 1, pp. 386–398, Jan. 2005.
](http://ieeexplore.ieee.org/document/1377521/)



