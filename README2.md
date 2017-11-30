# Introduction

LUT LDPC is a collection of software tools to design and test LDPC decoders based on discrete message passing decoding
using lookup tables (LUTs), referred to as LUT decoders, cf. [[1]](#Literature). It is mainly written in C++ and relies heavily on the [IT++](http://itpp.sourceforge.net/) for abstracting basic linear algebra and signal processing operations.
Consequently, the LUT decoders can easily be integrated into more complex communication systems including concatenated coding and/or modulation.

# Installation

## Requirements
The program has been tested on MacOSX and Linux, but with minor modification, will most likely also run on Windows. To compile you need static versions of the following libraries:
* IT++ and its dependencies
* boost
We favoured static over dynamic linking for portability, as we ran the compiled binaries on an inhomogeneous cluster of Linux hosts.
Dynamically linked versions could be obtained by changing the [Makefile](Makefile) appropriately.

To download and install the software open a terminal and enter
```
git clone --recursive https://github.com/mmeidlinger/lut_ldpc.git
cd lut_ldpc
make
```
This will compile and link the program as well as the patched IT++ library.
Note that all binaries and headers  are installed locally within the project directory.
Furthermore, you can pass the following options to `make`
* `BUILD_TYPE=Debug` (default: `Release`)


# Usage

LUT LDPC consists of several programs and scripts, whose usage we will discuss in what follows. In general, we want to emphasize that LUT LDPC is a research driven toolset. As such, its functionality should be used more like a library than via the included programs. Hence,
the ways to interact with the programs are limited and applying the software to specific and complex scenarios requires the user to write and compile their own programs. To assist users with that, LUT LDPC features [Doxygen](www.doxygen.org/) inline source code documentation.

## Designing LUT Decoders and testing Bit Error rate performance
### Running the simulation
`prog/ber_sim` is capable of designing decoders and conducting bit error rate (BER) Monte Carlo simulations.
```
; Sample parameter file for LDPC coded transmission over an BI-AWGN channel

; Some settings appear commented with their default values, unless the description of the setting
; explicitly states other default behavior. If you want to change
; them uncomment and edit them.

; Simulation settings
[Sim]
;  SNRs to simulate
   SNRdB    = 0:.5:4
;  Settings concerning the length and accuracy of a simulation
;  At most Nframes codewords are transmitted, however, transmission for a particular SNR
;  is terminated if the number exceeds Nfers. Furthermore, processing of additional SNR
;  points is terminated if either
;     * The bit error rate is below ber_min
;     * The frame error rate is below fer_min
   Nframes  = 1e2
;   Nfers    = 20
;   ber_min  = 1e-7
;   fer_min  = 1e-5
;  Random Seed Offset
;     In general the random seed is passed to the program when calling it. Here,
;     you can specify an offset
;   rand_seed_offset = 0
;  Custom name used to indentify the simulation
;   custom_name = 
;  Prefix for results
;   results_prefix = RES 
;  Settings regarding the search paths of the program. Keep in mind that the directories specified
;  here have to exist and the program needs to have permission to write to them. You can specify
;  them either relativ to the main directory of the program, but absolute paths also work fine.
;  As the search path for the parameter files has to be known to the program before loading param files,
;  it is specified by the BER_SIM_PARAMS_DIR preprocessor option. The default value ('.') can be
;  overridden at compile time.
;   results_dir = results
;   trees_dir = trees
;   codes_dir = codes
;  The codec is saved, if the random seed matches the value of save_codec. This is to avoid race conditions if multiple
;  instances of the program try to write the codec at the same time. Set save_codec to a negative number to avoid saving at all
;   save_codec=0
;  Codec filename and type. If these two setting are non-empty, the Codec is loaded from a file and all settings within the 
;  [LDPC], [BP] and [LUT] sections are ignored. This has to be supplied as an absolute path or relative to the working directory
;   codec_type = none ;BP, LUT
;   codec_filename = 


; Settings for the LDPC code. Further variables (codeword lendth, etc... 
; are determined upon loading the code)
[LDPC]
   parity_filename = rate0.50_dv02-17_dc08-09_lut_q4_N500
;  Wether a Generator matrix is used to encode random datawords or the all zero word is sent
   zero_codeword   = false 
;  Wether the Generator Matrix and the corresponding (possibly permuted) Parity Check Matrix that are generated
;  if LDPC.zero_codeword = true are saved to disk, potentially overwriting the original Parity check .alist file.
;   save_permuted = false
;  Wether a parity check is performed after each decoding iteration. Setting this to true is beneficial in case of a large number of iterattions.
;   parity_check_iter = true
   
; Settings for the BP Decoder. 
;[BP]
;  Maximum number of decoding operations. Warning: At least one option of this section
;  must be set for the program to recognize what kind of decoder to use.
;   max_iter = 30
;  Set the resolution of llrs
;     http://itpp.sourceforge.net/4.3.1/classitpp_1_1LLR__calc__unit.html#a531b0a4eca593e439a39ba9de42ee68c
;  Size of jac-log table. Set to 0 for min-sum decoding   
;   qllr_table_size = 300
;  Total number of bits for the qllr scaling factor determining the number of fractional bits
;   qllr_scale_res  = 12
;  Determines the spacing in between quantized LLR values
;   qllr_spacing_res = 7
;  Determines the internal resolution of the discrete LLRs. Defaults to 8*sizeof(QLLR)-4 if left empty
;   qllr_total_res  =


; Settings for a LUT based decoder. Warning: At least one option of this section
; must be set for the program to recognize what kind of decoder to use.
[LUT]
   max_iter = 50
;  Design SNR. Either this, or design_thr must be specified
;   design_SNRdB =
;  Design Noise threshold of an BIAWGN channel. If not specified, design_thr = sqrt(pow10(-design_SNRdB/10)/(2*LDPC_Code_Rate))
   design_thr = 0.88
;  Wether the check node update is LUT based or based on the min-sum rule
;   min_lut = true;
;  number of bits for channel LLR quantization
   qbits_channel = 4
;  number of bits per message. This assumes the resolution is uniform over iterations
   qbits_message_uniform = 4 
;  number of bits per message for each iteration. This has to be a vector of length max_iter
;  The default is to set it equal to qbits_message_uniform for all iteration. If you specify
;  the option qbits_message_uniform is overridden and does not have any impact
;  Example of non uniform resolution:
;   qbits_messages = 3 3 3 3 2 2 2 1
;  Set this vector to reuse Lookup-tables. The vector has to be 1 at position ii if at iteration ii, 
;  the LUT of iteration ii-1 is reused
;   reuse_lut = 0 0 1 0 0 1 0 0
;  Which tree structures to use for simulation options are
;     1. auto_bin_balanced: Automatically generate binary balanced trees
;     2. file: Load the tree structure frome a file within the 'trees' directory
;     3. auto_bin_high: Automatically generate a binary tree of maximum height
;     4. root_only: All leaf nodes are attached directly to the root node
;   tree_mode = auto_bin_balanced
;  Where to look for tree files
;   trees_dir = trees 
;   trees_filename =

```


### Evaluating the results

## Density evolution

## Generating codes



# Repository overview

### bin
Binary executables  corresponding to the sources in `prog`. Generated by `make` if need be; not under version control.

### build
Location of compiled object files  corresponding to the sources in `src`. Generated by `make` if need be; not under version control.

### codes
This is where we save LDPC parity check matrices, generator matrices and binary decoder objects, where only the parity check matrices (`.alist` files) should be under version control.

### doc
This is where the Doxygen documentation is generated

### ensembles
This is where LDPC ensembles (text files containing degree distributions are saved, cf. [ensembles/README.md](ensembles/README.md)

### include
Header directory

### itpp 
Location of the it++ library source, handled by a Git submodule, c.f. [IT++ Submodule](#it-submodule)

### lib
Directory for third party libraries (mainly itpp). Generated by `make` if need be; not under version control.

### params
Contains example parameter files for the program(s)

### prog
Contains the source code for executables with a `main()` function

### scripts
Various scripts, e.g. for evaluating and plotting BER simulation results or running simulations on an HTCondor cluster

### src
Contains the source code of the application modules, excluding executables with a `main()` function

### trees
This is where we save tree structures for the LUT based decoders



# Installation





# Contributing

## Changing the IT++ Submodul
Navigate into `cpp/itpp` ans issue `git status`. You can see that git complains about being in a detached head state. 
The reason for this is the way submodules work: The parent repository only saves a reference (commit-id) to the current commit of the
submodule. If you checkout the submodule, its head is set to this very commit. In order to create a local branch, type
``` git checkout -b ldpc-lut --track origin/ldpc-lut ```
This creates the local branch `ldpc-lut`, makes it point to the commit-id of the submodule and tells git that the
local branch `ldpc-lut` is supposed to track the remote branch `origin/ldpc-lut`.
You can now make changes to the submodule and push them upstream, given that you have write access to the submodule repo.
To get write access to the submodule containing the altered version of IT++ contact [Michael Meidlinger](mailto:michael@meidlinger.info).

# Referencing
If you use this software for your academic research, please consider referencing [our original contributions](#Literature)
# Literature
[[1] M. Meidlinger, A. Balatsoukas-Stimming, A. Burg, and G. Matz, “Quantized message passing for LDPC codes,” in Proc. 49th Asilomar Conf. Signals, Systems and Computers, Pacific Grove, CA, USA, Nov. 2015.](http://ieeexplore.ieee.org/document/7421419/)



