#!/bin/bash

module avail tool/idl
module load tool/idl-6.4
ulimit -s 6000000
ulimit -v unlimited

res='C90'
group='anthropogenic'
spec='CO'

for year in 1850
do
  cd /discover/nobackup/dgueyffi/modelE/aux

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  ./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_${spec}_anthropogenic_${year}_C90_Dec_2009.nc

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  ./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_${spec}_ships_${year}_C90_Dec_2009.nc

  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
  rm -f AR5-1.bat
  ./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 awb dom ind tra wst 
  echo ".com convert_${spec}.pro" >> ./AR5-1.bat
  echo ".run convert_${spec}.pro" >> ./AR5-1.bat
#now run the idl batch file:
  echo "exit" >> ./AR5-1.bat
  idl ./AR5-1.bat

  rm -f AR5-1.bat
  ./make_AR5_ships_${res}.ksh ${spec} ${year} ships C90_Dec_2009 shp
  echo ".com convert_${spec}.pro" >> ./AR5-1.bat
  echo ".run convert_${spec}.pro" >> ./AR5-1.bat
#now run the idl batch file:
  echo "exit" >> ./AR5-1.bat
  idl ./AR5-1.bat
done

for year in 1860 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990
do
  cd /discover/nobackup/dgueyffi/modelE/aux

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  ./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_${spec}_anthropogenic_${year}_C90_Dec_2009.nc

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  ./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_${spec}_ships_${year}_C90_Dec_2009.nc

  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
  rm -f AR5-1.bat
  ./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 awb dom ind tra wst ene
  echo ".com convert_${spec}.pro" >> ./AR5-1.bat
  echo ".run convert_${spec}.pro" >> ./AR5-1.bat
#now run the idl batch file:
  echo "exit" >> ./AR5-1.bat
  idl ./AR5-1.bat

  rm -f AR5-1.bat
  ./make_AR5_ships_${res}.ksh ${spec} ${year} ships C90_Dec_2009 shp
  echo ".com convert_${spec}.pro" >> ./AR5-1.bat
  echo ".run convert_${spec}.pro" >> ./AR5-1.bat
#now run the idl batch file:
  echo "exit" >> ./AR5-1.bat
  idl ./AR5-1.bat
done

for year in 2000
do
  cd /discover/nobackup/dgueyffi/modelE/aux

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  ./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_${spec}_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_${spec}_anthropogenic_${year}_C90_Dec_2009.nc

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  ./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_${spec}_ships_${year}_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_${spec}_ships_${year}_C90_Dec_2009.nc

  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
  rm -f AR5-1.bat
  ./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 agr awb dom ind tra wst ene slv
  echo ".com convert_${spec}.pro" >> ./AR5-1.bat
  echo ".run convert_${spec}.pro" >> ./AR5-1.bat
#now run the idl batch file:
  echo "exit" >> ./AR5-1.bat
  idl ./AR5-1.bat

  rm -f AR5-1.bat
  ./make_AR5_ships_${res}.ksh ${spec} ${year} ships C90_Dec_2009 shp
  echo ".com convert_${spec}.pro" >> ./AR5-1.bat
  echo ".run convert_${spec}.pro" >> ./AR5-1.bat
#now run the idl batch file:
  echo "exit" >> ./AR5-1.bat
  idl ./AR5-1.bat
done


for year in 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  cd /discover/nobackup/dgueyffi/modelE/aux

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_${spec}_decadalmonthlymean${year}_v1.nc.gz
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_${spec}_decadalmonthlymean${year}_v1.nc.gz .
  gunzip -f IPCC_GriddedBiomassBurningEmissions_${spec}_decadalmonthlymean${year}_v1.nc.gz

./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_${spec}_decadalmonthlymean${year}_v1.nc -out IPCC_GriddedBiomassBurningEmissions_${spec}_decadalmonthlymean${year}_C90_Dec_2009.nc

  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
rm -f AR5.bat
./make_AR5_program_BBURN_${res}_km.ksh ${spec} ${year} BiomassBurning C90_Dec_2009 grassfire forestfire
echo ".com convert_${spec}.pro" >> ./AR5.bat
echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
echo "exit" >> ./AR5.bat
idl ./AR5.bat

done

for year in 1850 1860 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
for src in awb dom ind tra wst shp
do
 fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1850-2000_${res} 1
done
done

for src in ene
do
fcop ./out_auto_C90/zero_annual_${res} ${spec}_${src}_AR5_1850-2000_${res}
for year in 1860 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1850-2000_${res} 1
done
done


for src in agr slv
do
fcop ./out_auto_C90/zero_annual_${res} ${spec}_${src}_AR5_1990-2000_${res}
for year in 2000
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1990-2000_${res} 1
done
done


for src in forestfire grassfire
do
for year in 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1890-2000_${res} 1
done
done
