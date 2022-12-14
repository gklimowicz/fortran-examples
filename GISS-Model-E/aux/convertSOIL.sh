#!/bin/bash
#
#
# This script performs the conversion to giss binary format for 
# fractions of silt, sand and clay for 2 layers: 0-30cm (tp file),
# 30-100cm (sb file)
#
# original data: binary file in non-compressed IDRISI format (.IMG) 
# The file contains bytes (from 0 to 255) in column major ordering.
# content: percentage of a given texture (e.g. silt)
# from:
# http://www.ngdc.noaa.gov/ecosys/cdroms/reynolds/reynolds/reynolds.htm
# C.A. Reynolds, T. J. Jackson, and W.J. Rawls. 1999. 
# Estimated Available Water Content from the FAO Soil Map of 
# the World, Global Soil Profile Databases, and Pedo-transfer Functions

# lanscape class
od -A n -t d1 -w9331200 landclss.img > landclss.txt
cat landclss.txt | awk 'BEGIN{n=1}{while(substr($0,n,5)){print substr($0,n,5);n+=5}}' > input
./bytetogiss
mv output.giss landclss.giss
 
#silt
cp /discover/nobackup/projects/giss/prod_input_files/silt_sb1.img .
od -A n -t d1 -w9331200 silt_sb1.img > silt_sb1.txt
cat silt_sb1.txt | awk 'BEGIN{n=1}{while(substr($0,n,5)){print substr($0,n,5);n+=5}}' > input
./bytetogiss
mv output.giss silt_sb1.giss

cp /discover/nobackup/projects/giss/prod_input_files/silt_tp1.img .
od -A n -t d1 -w9331200 silt_tp1.img > silt_tp1.txt
cat silt_tp1.txt | awk 'BEGIN{n=1}{while(substr($0,n,5)){print substr($0,n,5);n+=5}}' > input
./bytetogiss
mv output.giss silt_tp1.giss


# sand
cp /discover/nobackup/projects/giss/prod_input_files/sand_sb1.img .
od -A n -t d1 -w9331200 sand_sb1.img > sand_sb1.txt
cat sand_sb1.txt | awk 'BEGIN{n=1}{while(substr($0,n,5)){print substr($0,n,5);n+=5}}' > input
./bytetogiss
mv output.giss sand_sb1.giss

cp /discover/nobackup/projects/giss/prod_input_files/sand_tp1.img .
od -A n -t d1 -w9331200 sand_tp1.img > sand_tp1.txt
cat sand_tp1.txt | awk 'BEGIN{n=1}{while(substr($0,n,5)){print substr($0,n,5);n+=5}}' > input
./bytetogiss
mv output.giss sand_tp1.giss


# clay
cp /discover/nobackup/projects/giss/prod_input_files/clay_sb1.img .
od -A n -t d1 -w9331200 clay_sb1.img > clay_sb1.txt
cat clay_sb1.txt | awk 'BEGIN{n=1}{while(substr($0,n,5)){print substr($0,n,5);n+=5}}' > input
./bytetogiss
mv output.giss clay_sb1.giss

cp /discover/nobackup/projects/giss/prod_input_files/clay_tp1.img .
od -A n -t d1 -w9331200 clay_tp1.img > clay_tp1.txt
cat clay_tp1.txt | awk 'BEGIN{n=1}{while(substr($0,n,5)){print substr($0,n,5);n+=5}}' > input
./bytetogiss
mv output.giss clay_tp1.giss

ifort -convert big_endian convSOIL.f -o convSOIL
./convSOIL
