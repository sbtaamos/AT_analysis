#!/bin/bash

###############################################################################
# Usage: ./preprocess_AT.sh
###############################################################################

#source ~/.bashrc

# Amend for number of replicates and directories below:

for i in {1..5}; do

cd sim_$i/;

# Concatendate .xtc files for analysis:

gmx_avx trjcat -f md.*.part*.xtc -o md.PIP5K1A_popc_pops_pi4p_full_$i.xtc -nice -18;

gmx_avx trjconv -f md.PIP5K1A_popc_pops_pi4p_full_$i.xtc -dt 100 -o md.full_filtered_$i.xtc -nice 18

# Write a script that converts .xvg files to .dat files:

echo 'grep -v \# $1 | grep -v \@ > $2' >> xvg2dat.sh;

cd ../

printf "Sim $i concatenated."

done
