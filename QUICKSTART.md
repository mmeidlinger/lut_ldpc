# Quickstart Guide

* [Install dependencies for LUT-LDPC](#install-dependencies-for-lut-ldpc)
* [Download and Install LUT-LDPC and LUT-LDPC-VHDL](#download-and-install-lut-ldpc-and-lut-ldpc-vhdl)
* [Create preconfigured LDPC decoder object](#create-preconfigured-ldpc-decoder-object)
* [Create VHDL](#create-vhdl)
* [Install the free version of ModelSim](#install-the-free-version-of-modelsim)
   * [Download and install ModelSim free version](#download-and-install-modelsim-free-version)
   * [Install required libraries to run ModelSim](#install-required-libraries-to-run-modelsim)
* [Compile and run Testsuite](#compile-and-run-testsuite)

This file briefly goes through the shell commands for a complete [LUT-LDPC](README.md) decoder design flow
and the istallation of all required components.

We tested the following commands successfully on a fresh install of Ubuntu 17.10,
using only  software that is freely available.
So you can  go ahead and try out this tutorial on a Virtual Machime if you like.

## Install dependencies for LUT-LDPC
```shell
sudo apt-get install git build-essential cmake octave libboost-all-dev libfftw3-dev liblapack-dev libblas-dev
```

## Download and Install LUT-LDPC and LUT-LDPC-VHDL
```shell
cd ~/Desktop
git clone --recursive https://github.com/mmeidlinger/lut_ldpc.git
cd lut_ldpc
make
make install
```

## Create preconfigured LDPC decoder object
Design a LDPC decoder for the 10GBaseT code and run a short Bit Error Rate (BER) simulation, saving pairs of decoder inputs and outputs
into the text file `stimuli.txt`, cf. [here](https://github.com/mmeidlinger/lut_ldpc/blob/master/README.md#designing-lut-decoders-and-testing-bit-error-rate-performance) for more details. Those pairs can later be used to verify the functionality of
the VHDL decoder against the reference simulation.
```shell
bin/ber_sim -p params/ber.ini.regular.example > stimuli.txt
```
The resulting decoder object is saved under `results/RES_N2048_R0.841309_maxIter8_zcw0_frames20_minLUT/lut_codec.it`

## Create VHDL
Create VHDL code for an unrolled decoder based on the software decoder object produced by the above BER simulation:
```shell
cd lut_ldpc_vhdl/TopLevelDecoderGenerator
octave-cli --no-gui --verbose decoderGenerator.m
```

## Install the free version of ModelSim
ModelSim is used to verify the functionality of the generated VHDL code against the LUT-LDPC C++ simulation model.
To install ModelSim, we follow the instructions provided
[here](http://mattaw.blogspot.co.at/2014/05/making-modelsim-altera-starter-edition.html).

### Download and install ModelSim free version

```shell
cd /tmp
wget http://download.altera.com/akdlm/software/acdsinst/13.1/162/ib_installers/ModelSimSetup-13.1.0.162.run
chmod +x  ModelSimSetup-13.1.0.162.run
./ModelSimSetup-13.1.0.162.run
````
Install using the interactive installer. Select to install ModelSim to your home directory under `~/altera/13.1`.

### Install required libraries to run ModelSim
```shell
sudo dpkg --add-architecture i386
sudo apt-get update
sudo apt-get install gcc-multilib g++-multilib \
lib32z1 lib32stdc++6 lib32gcc1 \
expat:i386 fontconfig:i386 libfreetype6:i386 libexpat1:i386 libc6:i386 libgtk-3-0:i386 \
libcanberra0:i386 libpng-dev:i386 libice6:i386 libsm6:i386 libncurses5:i386 zlib1g:i386 \
libx11-6:i386 libxau6:i386 libxdmcp6:i386 libxext6:i386 libxft2:i386 libxrender1:i386 \
libxt6:i386 libxtst6:i386
```
```shell
cd /tmp
wget http://download.savannah.gnu.org/releases/freetype/freetype-2.4.12.tar.bz2
tar -xjvf freetype-2.4.12.tar.bz2
cd freetype-2.4.12
./configure --build=i686-pc-linux-gnu "CFLAGS=-m32" "CXXFLAGS=-m32" "LDFLAGS=-m32"
make -j4
mkdir -p ~/altera/13.1/modelsim_ase/lib32
cp objs/.libs/libfreetype.so* ~/altera/13.1/modelsim_ase/lib32
cd ~/altera/13.1/modelsim_ase
sed -i '/dir=`dirname $arg0`/a export LD_LIBRARY_PATH=${dir}/lib32' vco
ln -s linux linux_rh60
```

## Compile and run Testsuite

```shell
cd ~/Desktop/lut_ldpc/lut_ldpc_vhdl/ModelSim
./BIN_Decoder/compile.sh
./BIN_Decoder/run.sh
```
You can also open a ModelSim GUI with
```shell
~/altera/13.1/modelsim_ase/bin/vsim
```
