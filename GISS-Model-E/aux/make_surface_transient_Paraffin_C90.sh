#!/bin/bash

module avail tool/idl
module load tool/idl-6.4
ulimit -s 6000000
ulimit -v unlimited

res='C90'
group='anthropogenic'
spec='Paraffin'

for year in 1850 
do
  cd /discover/nobackup/dgueyffi/modelE/aux
  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz . 
  # gunzip -f IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 


  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz . 
  # gunzip -f IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst -w 0.03333333333333333333,0.01328021248339973439 IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst -w 1.0,1.0 propane+pentanes_anthropogenic.nc butanes+hexanes_anthropogenic.nc -o propane+pentanes+butanes+hexanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst -w 1.0,1.0 propane+pentanes+butanes+hexanes_anthropogenic.nc ethane+ketones_anthropogenic.nc -o paraffin_anthropogenic.nc


  ncflint -O -c -v emiss_shp -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_ships.nc
  ncflint -O -c -v emiss_shp -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_ships.nc
  ncflint -O -c -v emiss_shp -w 0.03333333333333333333,0.0 IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_ships.nc  #- note the weight == 0.0 there is no ketones_ships file
  ncflint -O -c -v emiss_shp -w 1.0,1.0 propane+pentanes_ships.nc butanes+hexanes_ships.nc -o propane+pentanes+butanes+hexanes_ships.nc
  ncflint -O -c -v emiss_shp -w 1.0,1.0 propane+pentanes+butanes+hexanes_ships.nc ethane+ketones_ships.nc -o paraffin_ships.nc

# then do regridding
./remap.pl -par ncregrid-ijl.par -in paraffin_ships.nc -out IPCC_emissions_Paraffin_ships_${year}_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in paraffin_anthropogenic.nc -out IPCC_emissions_Paraffin_anthropogenic_${year}_C90_Dec_2009.nc

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
  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz . 
  # gunzip -f IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 


  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz . 
  # gunzip -f IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra -w 0.03333333333333333333,0.01328021248339973439 IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra -w 1.0,1.0 propane+pentanes_anthropogenic.nc butanes+hexanes_anthropogenic.nc -o propane+pentanes+butanes+hexanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra -w 1.0,1.0 propane+pentanes+butanes+hexanes_anthropogenic.nc ethane+ketones_anthropogenic.nc -o paraffin_anthropogenic.nc


  ncflint -O -c -v emiss_shp -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_ships.nc
  ncflint -O -c -v emiss_shp -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_ships.nc
  ncflint -O -c -v emiss_shp -w 0.03333333333333333333,0.0 IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_ships.nc  #- note the weight == 0.0 there is no ketones_ships file
  ncflint -O -c -v emiss_shp -w 1.0,1.0 propane+pentanes_ships.nc butanes+hexanes_ships.nc -o propane+pentanes+butanes+hexanes_ships.nc
  ncflint -O -c -v emiss_shp -w 1.0,1.0 propane+pentanes+butanes+hexanes_ships.nc ethane+ketones_ships.nc -o paraffin_ships.nc


# then do regridding
./remap.pl -par ncregrid-ijl.par -in paraffin_ships.nc -out IPCC_emissions_Paraffin_ships_${year}_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in paraffin_anthropogenic.nc -out IPCC_emissions_Paraffin_anthropogenic_${year}_C90_Dec_2009.nc

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
  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz . 
  # gunzip -f IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 


  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz . 
  # gunzip -f IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 


  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra,emiss_awb -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra,emiss_awb -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra,emiss_awb -w 0.03333333333333333333,0.01328021248339973439 IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra,emiss_awb -w 1.0,1.0 propane+pentanes_anthropogenic.nc butanes+hexanes_anthropogenic.nc -o propane+pentanes+butanes+hexanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra,emiss_awb -w 1.0,1.0 propane+pentanes+butanes+hexanes_anthropogenic.nc ethane+ketones_anthropogenic.nc -o paraffin_anthropogenic.nc
  

ncflint -O -c -v emiss_shp -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_ships.nc
ncflint -O -c -v emiss_shp -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_ships.nc
ncflint -O -c -v emiss_shp -w 0.03333333333333333333,0.0 IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_ships.nc  #- note the weight == 0.0 there is no ketones_ships file
ncflint -O -c -v emiss_shp -w 1.0,1.0 propane+pentanes_ships.nc butanes+hexanes_ships.nc -o propane+pentanes+butanes+hexanes_ships.nc
ncflint -O -c -v emiss_shp -w 1.0,1.0 propane+pentanes+butanes+hexanes_ships.nc ethane+ketones_ships.nc -o paraffin_ships.nc

# then do regridding
./remap.pl -par ncregrid-ijl.par -in paraffin_ships.nc -out IPCC_emissions_Paraffin_ships_${year}_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in paraffin_anthropogenic.nc -out IPCC_emissions_Paraffin_anthropogenic_${year}_C90_Dec_2009.nc

cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
rm -f AR5.bat
./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 dom ind wst ene tra awb
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

for year in 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  cd /discover/nobackup/dgueyffi/modelE/aux

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean${year}_v1.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean${year}_v1.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean${year}_v1.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean${year}_v1.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean${year}_v1.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean${year}_v1.nc.gz

ncflint -O -c -v grassfire,forestfire -w 0.02272727272727272727,0.01388888888888888888 IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean${year}_v1.nc IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean${year}_v1.nc -o propane+pentanes_Biomass.nc
ncflint -O -c -v grassfire,forestfire -w 0.01730103806228373702,0.00936329588014981273 IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean${year}_v1.nc IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean${year}_v1.nc -o butanes+hexanes_Biomass.nc
ncflint -O -c -v grassfire,forestfire -w 0.03333333333333333333,0.01328021248339973439 IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean${year}_v1.nc IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean${year}_v1.nc -o ethane+ketones_Biomass.nc
ncflint -O -c -v grassfire,forestfire -w 1.0,1.0 propane+pentanes_Biomass.nc butanes+hexanes_Biomass.nc -o propane+pentanes+butanes+hexanes_Biomass.nc
ncflint -O -c -v grassfire,forestfire -w 1.0,1.0 propane+pentanes+butanes+hexanes_Biomass.nc ethane+ketones_Biomass.nc -o paraffin_Biomass.nc

# then do regridding
./remap.pl -par ncregrid-ijl.par -in paraffin_Biomass.nc -out IPCC_GriddedBiomassBurningEmissions_Paraffin_decadalmonthlymean${year}_C90_Dec_2009.nc

cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro

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
  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz . 
  # gunzip -f IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 


  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz . 
  # gunzip -f IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean${year}_v1.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean${year}_v1.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean${year}_v1.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean${year}_v1.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean${year}_v1.nc.gz

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean${year}_v1.nc.gz
  # cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean${year}_v1.nc.gz .
  # gunzip -f IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean${year}_v1.nc.gz


  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra,emiss_awb,emiss_agr -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra,emiss_awb,emiss_agr -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra,emiss_awb,emiss_agr -w 0.03333333333333333333,0.01328021248339973439 IPCC_emissions_ethane_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra,emiss_awb,emiss_agr -w 1.0,1.0 propane+pentanes_anthropogenic.nc butanes+hexanes_anthropogenic.nc -o propane+pentanes+butanes+hexanes_anthropogenic.nc
  ncflint -O -c -v emiss_dom,emiss_ind,emiss_wst,emiss_ene,emiss_tra,emiss_awb,emiss_agr -w 1.0,1.0 propane+pentanes+butanes+hexanes_anthropogenic.nc ethane+ketones_anthropogenic.nc -o paraffin_anthropogenic.nc
  
  ncflint -O -c -v emiss_shp -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_ships.nc
  ncflint -O -c -v emiss_shp -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_ships.nc
  ncflint -O -c -v emiss_shp -w 0.03333333333333333333,0.0 IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ethane_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_ships.nc  #- note the weight == 0.0 there is no ketones_ships file
  ncflint -O -c -v emiss_shp -w 1.0,1.0 propane+pentanes_ships.nc butanes+hexanes_ships.nc -o propane+pentanes+butanes+hexanes_ships.nc
  ncflint -O -c -v emiss_shp -w 1.0,1.0 propane+pentanes+butanes+hexanes_ships.nc ethane+ketones_ships.nc -o paraffin_ships.nc

ncflint -O -c -v grassfire,forestfire -w 0.02272727272727272727,0.01388888888888888888 IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean${year}_v1.nc IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean${year}_v1.nc -o propane+pentanes_Biomass.nc
ncflint -O -c -v grassfire,forestfire -w 0.01730103806228373702,0.00936329588014981273 IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean${year}_v1.nc IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean${year}_v1.nc -o butanes+hexanes_Biomass.nc
ncflint -O -c -v grassfire,forestfire -w 0.03333333333333333333,0.01328021248339973439 IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean${year}_v1.nc IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean${year}_v1.nc -o ethane+ketones_Biomass.nc
ncflint -O -c -v grassfire,forestfire -w 1.0,1.0 propane+pentanes_Biomass.nc butanes+hexanes_Biomass.nc -o propane+pentanes+butanes+hexanes_Biomass.nc
ncflint -O -c -v grassfire,forestfire -w 1.0,1.0 propane+pentanes+butanes+hexanes_Biomass.nc ethane+ketones_Biomass.nc -o paraffin_Biomass.nc

# then do regridding
./remap.pl -par ncregrid-ijl.par -in paraffin_anthropogenic.nc -out IPCC_emissions_Paraffin_anthropogenic_${year}_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in paraffin_ships.nc -out IPCC_emissions_Paraffin_ships_${year}_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in paraffin_Biomass.nc -out IPCC_GriddedBiomassBurningEmissions_Paraffin_decadalmonthlymean${year}_C90_Dec_2009.nc

cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
rm -f AR5.bat
./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 dom ind wst ene tra awb agr 
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



for year in 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  cd /discover/nobackup/dgueyffi/modelE/aux

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  # dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz  
  # cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  # gunzip -f IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc.gz 

  ncflint -O -c -v emiss_slv -w 0.00936329588014981273,0.01328021248339973439 IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ketones_anthropogenic_${year}_0.5x0.5_v1_20_04_2009.nc -o paraffin_slv.nc

# then do regridding
./remap.pl -par ncregrid-ijl.par -in paraffin_slv.nc -out IPCC_emissions_Paraffin_anthropogenic_${year}_C90_Dec_2009.nc

cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
rm -f AR5.bat
./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 slv
echo ".com convert_${spec}.pro" >> ./AR5.bat
echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
echo "exit" >> ./AR5.bat
idl ./AR5.bat

done 


for src in dom ind wst shp
do
for year in 1850 1860 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
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



for src in slv
do          
fcop ./out_auto_C90/zero_annual_${res} ${spec}_${src}_AR5_1860-2000_${res}
for year in 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1860-2000_${res} 1
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
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1900-2000_${res} 1
done
done
