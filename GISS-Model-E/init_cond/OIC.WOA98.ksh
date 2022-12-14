#!/bin/ksh  
USAGE="usage: OIC.WOA98.ksh" 
#
# This script complies OIC.WOA98.f, and then computes the ocean IC. 
#  The original data are in the dir. /discover/nobackup/projects/giss/OBS; 
#  so place the data in this dir. or change the path for OBS.
#  Right now the OIC.WOA98.f is set up for IM=144,JM=90, and for 32 layer 
#  data to create. If you need to create the data for different dimensions, 
#  then read carefully the comments in the beginning of OIC.WOA98.f and 
#  make all necessary changes.
#
export IFDIR=/discover/nobackup/projects/giss/prod_input_files/
export OBS=/discover/nobackup/projects/giss/OBS/
#
ifort -convert big_endian -o OIC.WOA98.exe OIC.WOA98.f
OIC.WOA98.exe
#
#mv OIC.WOA98.2HX2.L32.D1201 $IFDIR 
# 
rm OIC.WOA98.exe






