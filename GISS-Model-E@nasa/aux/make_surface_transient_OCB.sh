#!/bin/bash

module avail tool/idl
module load tool/idl-6.4
ulimit -s 6000000
ulimit -v unlimited

res='C90'
spec='OCB'

for year in 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  cd /discover/nobackup/dgueyffi/modelE/aux

  dmget /archive/g08/dmkoch/ftp/AR5/IPCC_GriddedBiomassBurningEmissions_${spec}_decadalmonthlymean${year}_v1.nc
  cp /archive/g08/dmkoch/ftp/AR5/IPCC_GriddedBiomassBurningEmissions_${spec}_decadalmonthlymean${year}_v1.nc .

 ./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_${spec}_decadalmonthlymean${year}_v1.nc -out IPCC_GriddedBiomassBurningEmissions_${spec}_decadalmonthlymean${year}_C90_Dec_2009.nc

  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
rm -f AR5-OCB.dat
./make_AR5_program_BBURN_${res}.ksh ${spec} ${year} BiomassBurning C90_Dec_2009 grassfire forestfire
echo ".com convert_${spec}.pro" >> ./AR5-OCB.dat
echo ".run convert_${spec}.pro" >> ./AR5-OCB.dat
#now run the idl batch file:
echo "exit" >> ./AR5-OCB.dat
idl ./AR5-OCB.dat

done

for src in forestfire grassfire
do
for year in 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1890-2000_${res} 1
done
done
