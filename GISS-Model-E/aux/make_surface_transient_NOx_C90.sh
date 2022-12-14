#!/bin/ksh

res='C90'
group='anthropogenic'
spec='NOx'

for year in 1850 1860 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  cd /discover/nobackup/dgueyffi/modelE/aux
  rm -f NOx_anthropo.nc
  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_NO_anthropogenic_${year}_0.5x0.5_v1_07_05_2009.nc.gz
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_NO_anthropogenic_${year}_0.5x0.5_v1_07_05_2009.nc.gz .
  gunzip -f IPCC_emissions_NO_anthropogenic_${year}_0.5x0.5_v1_07_05_2009.nc.gz

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_NO_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz
  cp  /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_NO_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz .
  gunzip -f IPCC_emissions_NO_ships_${year}_0.5x0.5_v1_20_04_2009.nc.gz

  ncflint -O -c -v emiss_agr,emiss_awb,emiss_dom,emiss_ene,emiss_ind,emiss_tra,emiss_wst -w .46666666666666666666,0.0 IPCC_emissions_NO_anthropogenic_${year}_0.5x0.5_v1_07_05_2009.nc IPCC_emissions_NO_anthropogenic_${year}_0.5x0.5_v1_07_05_2009.nc -o NOx_anthropo.nc
  
  ncflint -O -c -v emiss_shp -w .46666666666666666666,0.0 IPCC_emissions_NO_ships_${year}_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_NO_ships_${year}_0.5x0.5_v1_20_04_2009.nc -o NOx_ships.nc

  ./remap.pl -par ncregrid-ijl.par -in NOx_anthropo.nc -out IPCC_emissions_${spec}_anthropogenic_${year}_C90_Dec_2009.nc
  ./remap.pl -par ncregrid-ijl.par -in NOx_ships.nc -out IPCC_emissions_${spec}_ships_${year}_C90_Dec_2009.nc

  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
  rm -f AR5.bat
  ./make_AR5_program_${res}.ksh ${spec} ${year} ${group} C90_Dec_2009 ene dom ind wst agr awb tra
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

  rm -f AR5.bat
  ./make_AR5_ships_${res}.ksh ${spec} ${year} ships C90_Dec_2009 shp
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat

done

for year in 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  cd /discover/nobackup/dgueyffi/modelE/aux

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_NOx_decadalmonthlymean${year}_v1.nc.gz
  cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_GriddedBiomassBurningEmissions_NOx_decadalmonthlymean${year}_v1.nc.gz .
  gunzip -f IPCC_GriddedBiomassBurningEmissions_NOx_decadalmonthlymean${year}_v1.nc.gz

  ncflint -O -c -v forestfire,grassfire -w .46666666666666666666,0.0 IPCC_GriddedBiomassBurningEmissions_NOx_decadalmonthlymean${year}_v1.nc IPCC_GriddedBiomassBurningEmissions_NOx_decadalmonthlymean${year}_v1.nc -o NOx_Biomass.nc
./remap.pl -par ncregrid-ijl.par -in NOx_Biomass.nc -out IPCC_GriddedBiomassBurningEmissions_${spec}_decadalmonthlymean${year}_C90_Dec_2009.nc

  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
rm -f AR5.bat
./make_AR5_program_BBURN_${res}_km.ksh ${spec} ${year} BiomassBurning C90_Dec_2009 grassfire forestfire
echo ".com convert_${spec}.pro" >> ./AR5.bat
echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
echo "exit" >> ./AR5.bat
idl ./AR5.bat

done

for year in 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do 
  cd /discover/nobackup/dgueyffi/modelE/aux

  dmget -v /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_NO_aircraft_${year}_0.5x0.5_v1_23_07_2009.nc.gz
  cp /nfs3m/archive/g08/gfaluveg/AR5_emissions/RAW/${year}/IPCC_emissions_NO_aircraft_${year}_0.5x0.5_v1_23_07_2009.nc.gz .
  gunzip -f IPCC_emissions_NO_aircraft_${year}_0.5x0.5_v1_23_07_2009.nc.gz

  ncflint -O -c -v emiss_air -w .46666666666666666666,0.0 IPCC_emissions_NO_aircraft_${year}_0.5x0.5_v1_23_07_2009.nc IPCC_emissions_NO_aircraft_${year}_0.5x0.5_v1_23_07_2009.nc -o NOx_air.nc
  
  ./remap.pl -par ncregrid-ijkl.par -in NOx_air.nc -out IPCC_emissions_${spec}_aircraft_${year}_C90_Dec_2009.nc

  cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
  rm -f AR5.bat
  ./make_AR5_program_aircraft_${res}.ksh ${spec} ${year} aircraft C90_Dec_2009 air
  echo ".com convert_${spec}.pro" >> ./AR5.bat
  echo ".run convert_${spec}.pro" >> ./AR5.bat
#now run the idl batch file:
  echo "exit" >> ./AR5.bat
  idl ./AR5.bat
done

for year in 1850 1860 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
for src in agr awb dom ene ind tra wst shp 
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1850-2000_${res} 1
done
done

for year in 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
for src in air
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1910-2000_${res}
done
done

for src in forestfire grassfire
do
for year in 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
do
  fcop ./out_auto_C90/${year}/${res}/${spec}_${src}_AR5_${year}_${res}_h ${spec}_${src}_AR5_1890-2000_${res} 1
done
done


