#!/bin/bash

module avail tool/idl
module load tool/idl-6.4
ulimit -s 6000000
ulimit -v unlimited

res='C90'
group='anthropogenic'
spec='Alkenes'

for year in 1850
do
  cd /discover/nobackup/dgueyffi/modelE/aux
  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst -w 1.0,0.03571428571428571428 propene+oalkenes_anthropogenic.nc IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o alkenes_anthropogenic.nc
  
  ncflint -O -c -v emiss_shp -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_ships.nc 
  ncflint -O -c -v emiss_shp -w 1.0,0.03571428571428571428 propene+oalkenes_ships.nc IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o alkenes_ships.nc 

./remap.pl -par ncregrid-ijl.par -in alkenes_anthropogenic.nc -out IPCC_emissions_Alkenes_anthropogenic_${year}_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in alkenes_ships.nc -out IPCC_emissions_Alkenes_ships_${year}_C90_Dec_2009.nc

  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
  rm -f AR5.bat
  ./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 dom ind wst
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

  rm -f AR5.bat
  ./make_AR5_ships_${res}_km.ksh ${spec} ${year} ships C90_Dec_2009 shp
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

done


for year in 1860 1870 1880 1890
do
  cd /discover/nobackup/dgueyffi/modelE/aux
  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  ncflint -O -c -v emiss_dom,emiss_ene,emiss_ind,emiss_tra,emiss_wst -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ene,emiss_ind,emiss_tra,emiss_wst -w 1.0,0.03571428571428571428 propene+oalkenes_anthropogenic.nc IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o alkenes_anthropogenic.nc
  
  ncflint -O -c -v emiss_shp -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_ships.nc 
  ncflint -O -c -v emiss_shp -w 1.0,0.03571428571428571428 propene+oalkenes_ships.nc IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o alkenes_ships.nc 

  ./remap.pl -par ncregrid-ijl.par -in alkenes_anthropogenic.nc -out IPCC_emissions_Alkenes_anthropogenic_${year}_C90_Dec_2009.nc
  ./remap.pl -par ncregrid-ijl.par -in alkenes_ships.nc -out IPCC_emissions_Alkenes_ships_${year}_C90_Dec_2009.nc
  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
  rm -f AR5.bat
  ./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 dom ind wst ene tra 
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

  rm -f AR5.bat
  ./make_AR5_ships_${res}_km.ksh ${spec} ${year} ships C90_Dec_2009 shp
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

done


for year in 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990
do
  cd /discover/nobackup/dgueyffi/modelE/aux
  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean${year}_v1.nc.gz 
  cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean${year}_v1.nc.gz .
  gunzip -f IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean${year}_v1.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean${year}_v1.nc.gz
  cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean${year}_v1.nc.gz .
  gunzip -f IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean${year}_v1.nc.gz 

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean${year}_v1.nc.gz
  cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean${year}_v1.nc.gz .
  gunzip -f IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean${year}_v1.nc.gz


  ncflint -O -c -v emiss_dom,emiss_ene,emiss_ind,emiss_awb,emiss_tra,emiss_wst -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ene,emiss_ind,emiss_awb,emiss_tra,emiss_wst -w 1.0,0.03571428571428571428 propene+oalkenes_anthropogenic.nc IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o alkenes_anthropogenic.nc
  
  ncflint -O -c -v emiss_shp -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_ships.nc 
  ncflint -O -c -v emiss_shp -w 1.0,0.03571428571428571428 propene+oalkenes_ships.nc IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o alkenes_ships.nc 

  ncflint -O -c -v grassfire,forestfire -w 0.02380952380952380952,0.01492537313432835820 IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean${year}_v1.nc IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean${year}_v1.nc -o propene+oalkenes_Biomass.nc 
  ncflint -O -c -v grassfire,forestfire -w 1.0,0.03571428571428571428 propene+oalkenes_Biomass.nc IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean${year}_v1.nc -o alkenes_Biomass.nc

./remap.pl -par ncregrid-ijl.par -in alkenes_anthropogenic.nc -out IPCC_emissions_Alkenes_anthropogenic_${year}_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in alkenes_ships.nc -out IPCC_emissions_Alkenes_ships_${year}_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in alkenes_Biomass.nc -out IPCC_GriddedBiomassBurningEmissions_Alkenes_decadalmonthlymean${year}_C90_Dec_2009.nc
  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro

  rm -f AR5.bat
  ./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 dom ind wst awb ene tra
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

  rm -f AR5.bat
  ./make_AR5_ships_${res}_km.ksh ${spec} ${year} ships C90_Dec_2009 shp
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

  rm -f AR5.bat
  ./make_AR5_program_BBURN_${res}_km.ksh ${spec} ${year} BiomassBurning C90_Dec_2009 grassfire forestfire
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

done


for year in 2000
do
  cd /discover/nobackup/dgueyffi/modelE/aux
  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean${year}_v1.nc.gz 
  cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean${year}_v1.nc.gz .
  gunzip -f IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean${year}_v1.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean${year}_v1.nc.gz
  cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean${year}_v1.nc.gz .
  gunzip -f IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean${year}_v1.nc.gz 

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean${year}_v1.nc.gz
  cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean${year}_v1.nc.gz .
  gunzip -f IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean${year}_v1.nc.gz

  ncflint -O -c -v emiss_dom,emiss_ene,emiss_ind,emiss_awb,emiss_agr,emiss_tra,emiss_wst -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ene,emiss_ind,emiss_awb,emiss_agr,emiss_tra,emiss_wst -w 1.0,0.03571428571428571428 propene+oalkenes_anthropogenic.nc IPCC_emissions_ethene_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o alkenes_anthropogenic.nc

  ncflint -O -c -v emiss_shp -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_ships.nc 
  ncflint -O -c -v emiss_shp -w 1.0,0.03571428571428571428 propene+oalkenes_ships.nc IPCC_emissions_ethene_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o alkenes_ships.nc 

ncflint -O -c -v grassfire,forestfire -w 0.02380952380952380952,0.01492537313432835820 IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean${year}_v1.nc IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean${year}_v1.nc -o propene+oalkenes_Biomass.nc 
ncflint -O -c -v grassfire,forestfire -w 1.0,0.03571428571428571428 propene+oalkenes_Biomass.nc IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean${year}_v1.nc -o alkenes_Biomass.nc

./remap.pl -par ncregrid-ijl.par -in alkenes_anthropogenic.nc -out IPCC_emissions_Alkenes_anthropogenic_${year}_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in alkenes_ships.nc -out IPCC_emissions_Alkenes_ships_${year}_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in alkenes_Biomass.nc -out IPCC_GriddedBiomassBurningEmissions_Alkenes_decadalmonthlymean${year}_C90_Dec_2009.nc
  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
  rm -f AR5.bat
  ./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 dom ind wst awb agr ene tra
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

  rm -f AR5.bat
  ./make_AR5_ships_${res}_km.ksh ${spec} ${year} ships C90_Dec_2009 shp
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

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
for src in dom ind wst shp
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1850-2000_${res} 1
done
done

for src in ene tra
do
fcop ./out_auto_C90/zero_annual_${res} ${spec}_${src}_AR5_1850-2000_${res}
for year in 1860 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1850-2000_${res} 1
done
done

for src in awb 
do
fcop ./out_auto_C90/zero_annual_${res} ${spec}_${src}_AR5_1890-2000_${res}
for year in 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1890-2000_${res} 1
done
done

for src in agr
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