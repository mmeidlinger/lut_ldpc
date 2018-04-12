# Entire designflow to produce the VHDL for an unrolled LUT decoder. Tested on Ubuntu 17.10.1 on 04-12-2018 by Michael Meidlinger (michael@meidlinger.info)

# Install dependencies for LUT-LDPC
sudo apt-get install git build-essential cmake octave libboost-all-dev libfftw3-dev liblapack-dev libblas-dev
# Download and Install LUT-LDPC and LUT-LDPC-VHDL:
cd ~/Desktop
git clone --recursive https://github.com/mmeidlinger/lut_ldpc.git
cd lut_ldpc
make
make install
# Design a LDPC decoder for the 10GBaseT code and run a short Bit Error Rate (BER) simulation:
bin/ber_sim -p params/ber.ini.regular.example > stimuli.txt
# Create VHDL for an unrolled decoder based on the software decoder object produced by the above BER simulation:
cd lut_ldpc_vhdl/TopLevelDecoderGenerator
octave-cli --no-gui --verbose decoderGenerator.m 


# Install the free version of ModelSim to verify the VHDL code against the simulation model. We follow the instructions provided at
# http://mattaw.blogspot.co.at/2014/05/making-modelsim-altera-starter-edition.html

# Donload and install ModelSim free version
cd /tmp
wget http://download.altera.com/akdlm/software/acdsinst/13.1/162/ib_installers/ModelSimSetup-13.1.0.162.run
chmod +x  ModelSimSetup-13.1.0.162.run
./ModelSimSetup-13.1.0.162.run
# Install using the interactive installer. Select to install ModelSim to your home directory under ~/altera/13.1

# Install required libraries to run ModelSim
sudo dpkg --add-architecture i386
sudo apt-get update
sudo apt-get install gcc-multilib g++-multilib \
lib32z1 lib32stdc++6 lib32gcc1 \
expat:i386 fontconfig:i386 libfreetype6:i386 libexpat1:i386 libc6:i386 libgtk-3-0:i386 \
libcanberra0:i386 libpng-dev:i386 libice6:i386 libsm6:i386 libncurses5:i386 zlib1g:i386 \
libx11-6:i386 libxau6:i386 libxdmcp6:i386 libxext6:i386 libxft2:i386 libxrender1:i386 \
libxt6:i386 libxtst6:i386
wget http://download.savannah.gnu.org/releases/freetype/freetype-2.4.12.tar.bz2
tar -xjvf freetype-2.4.12.tar.bz2
cd freetype-2.4.12
./configure --build=i686-pc-linux-gnu "CFLAGS=-m32" "CXXFLAGS=-m32" "LDFLAGS=-m32"
make -j4	
mkdir -p ~/altera/13.1/modelsim_ase/lib32
cp objs/.libs/libfreetype.so* ~/altera/13.1/modelsim_ase/lib32
cd ~/altera/13.1/modelsim_ase
sed '/dir=`dirname $arg0`/a export LD_LIBRARY_PATH=${dir}/lib32' vco
ln -s linux linux_rh60

# Compile and run Testsuite
cd ~/Desktop/lut_ldpc/lut_ldpc_vhdl/ModelSim
./BIN_Decoder/compile.sh 
./BIN_Decoder/run.sh

# Open a ModelSim GUI
~/altera/13.1/modelsim_ase/bin/vsim

