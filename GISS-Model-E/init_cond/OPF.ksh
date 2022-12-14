#!/bin/ksh  
USAGE="usage: OPF.ksh" 
#
# This script complies OPF.f, and then computes the polar filter data.  
#  The input data file Z144X90N_nocasp.1 is taken from the hard-coded dir. 
#  IFDIR=/discover/nobackup/projects/giss/prod_input_files/ on Discover; 
#  so change the path for IFDIR.
#  Right now the OPF.f is set up for IM=144,JM=90, and 32 layer data 
#  data to create. If you need to create the data for different dimensions, 
#  then read carefully the comments in the beginning of OPF.f and 
#  make all necessary changes.
#
export IFDIR=/discover/nobackup/projects/giss/prod_input_files/

ifort -convert big_endian -o OPF.exe OPF.f
OPF.exe
#
#mv OPF.E2HX2.L32 $IFDIR 
# 
rm OPF.exe






