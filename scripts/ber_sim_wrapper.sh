#!/bin/sh

# The condor executable is usually transferred over the network to the machines executing the jobs. Using this script rather than a binary thus limits the communication  overhead.

params_file=$1
rand_seed=$2
base_dir=$3
custom_name=$4

# Call to the actual program.
$base_dir/bin/ber_sim -p $params_file -s $rand_seed -b $base_dir -c $custom_name
