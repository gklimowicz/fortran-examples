#!/bin/bash
# Regridding IPCC emissions files to the C90 cubed sphere grid
# You must run this script on a DALI node 
#cp /discover/nobackup/projects/giss/prod_input_files/IPCC_emissions*_1850_0.5x0.5_v1_20_04_2009.nc .

#- cleanup
rm -f paraffin_anthropogenic.nc propane+pentanes_anthropogenic.nc propane+pentanes_Biomass.nc propane+pentanes+butanes+hexanes_anthropogenic.nc
rm -f propane+pentanes+butanes+hexanes_ships.nc propane+pentanes_ships.nc propene+oalkenes_anthropogenic.nc propene+oalkenes_Biomass.nc
rm -f propene+oalkenes_ships.nc alkenes_anthropogenic.nc alkenes_Biomass.nc alkenes_ships.nc 
rm -f butanes+hexanes_* ethane+ketones_*nc propane+pentanes+butanes+hexanes_Biomass.nc paraffin_Biomass.nc
rm -f paraffin_ships.nc NOx*.nc

#-- vegetation density
./remap.pl -par regrida2x2.5.par -in veg_dense_2x2.5 -out veg_dense_C90_from_2x2.5

#-- codirect
./remap.pl -par ncregrid-ijl.par -in HTAP_surface_co_2001_aggregated.0.5x0.5.nc -out HTAP_surface_co_2001_aggregated_C90.nc 

#- some dust files
./remap.pl -par regridb2x2.5.par -in GriniZender_DustSources_144x90 -out GriniZender_DustSources_C90_from_144x90 
./remap.pl -par regridb2x2.5.par -in Tegen_DustSources_144x90 -out Tegen_DustSources_C90_from_144x90
#- Tegen from high res 0.5x0.5 source file
./remap.pl -par ncregrid-tegen.par -in Tegen_source_redresse.nc -out Tegen02_C90.nc

#-- interactive wetland
./remap.pl -par regrida2x2.5.par -in beta_p_ch4_4x5_2x2.5gf -out beta_p_ch4_4x5_2x2.5gf_C90
./remap.pl -par regrida2x2.5.par -in alpha_t_ch4_4x5_2x2.5gf -out alpha_t_ch4_4x5_2x2.5gf_C90
./remap.pl -par regrida2x2.5.par -in ncep_g1temp_2x2.5 -out ncep_g1temp_2x2.5_C90
./remap.pl -par regrida2x2.5.par -in ncep_prec_w_2wk_lag_2x2.5 -out ncep_prec_w_2wk_lag_2x2.5_C90

#-- interactive fire emissions
./remap.pl -par regrida2x2.5.par -in OC_emmissions_per_fire_count_2x2.5 -out OC_emmissions_per_fire_count_C90_from_2x2.5
./remap.pl -par regrida2x2.5.par -in BC_emmissions_per_fire_count_2x2.5 -out BC_emmissions_per_fire_count_C90_from_2x2.5
./remap.pl -par regrida2x2.5.par -in Paraffin_emmissions_per_fire_count_2x2.5 -out Paraffin_emmissions_per_fire_count_C90_from_2x2.5
./remap.pl -par regrida2x2.5.par -in Alkenes_emmissions_per_fire_count_2x2.5 -out Alkenes_emmissions_per_fire_count_C90_from_2x2.5
./remap.pl -par regrida2x2.5.par -in CH4_emmissions_per_fire_count_2x2.5 -out CH4_emmissions_per_fire_count_C90_from_2x2.5
./remap.pl -par regrida2x2.5.par -in CO_emmissions_per_fire_count_2x2.5 -out CO_emmissions_per_fire_count_C90_from_2x2.5
./remap.pl -par regrida2x2.5.par -in NOx_emmissions_per_fire_count_2x2.5 -out NOx_emmissions_per_fire_count_C90_from_2x2.5

#-- Gas tracers - initial conditions -
./remap.pl -par ICgastracers.par -in N2O_IC_M23_4x5_6.17_conc_2x2.5_conc -out N2O_IC_from_2x2.5_C90_Dec_2009 
./remap.pl -par ICgastracers.par -in CFC_IC_M23_4x5_6.17_conc_2x2.5_conc -out CFC_IC_from_2x2.5_C90_Dec_2009 
./remap.pl -par ICgastracers.par -in CH4_IC_M23_4x5_6.17_conc_2x2.5_conc -out CH4_IC_from_2x2.5_C90_Dec_2009 
./remap.pl -par ICgastracers.par -in Ox_init_cond_M23_4x5_conc_2x2.5_conc -out Ox_IC_from_2x2.5_C90_Dec_2009 
./remap.pl -par ICgastracers.par -in CO_init_cond_M23_conc_2x2.5_conc -out CO_IC_from_2x2.5_C90_Dec_2009 
./remap.pl -par ICOxref.par -in O3ref_O3JDAY_1850_182.dat -out O3ref_O3JDAY_1850_from_2x2.5_C90_Dec_2009  #check this one
./remap.pl -par ICsulfate.par -in sulfate_pi_fakeM23_M_SA_2x2.5gf -out sulfate_from_2x2.5_C90_Dec_2009

#-- Ozone
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1850 -out jan2010_o3_shindell_144x90x49x12_1850_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1870 -out jan2010_o3_shindell_144x90x49x12_1870_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1890 -out jan2010_o3_shindell_144x90x49x12_1890_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1910 -out jan2010_o3_shindell_144x90x49x12_1910_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1930 -out jan2010_o3_shindell_144x90x49x12_1930_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1940 -out jan2010_o3_shindell_144x90x49x12_1940_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1950 -out jan2010_o3_shindell_144x90x49x12_1950_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1960 -out jan2010_o3_shindell_144x90x49x12_1960_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1970 -out jan2010_o3_shindell_144x90x49x12_1970_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1980 -out jan2010_o3_shindell_144x90x49x12_1980_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_1990 -out jan2010_o3_shindell_144x90x49x12_1990_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_2000 -out jan2010_o3_shindell_144x90x49x12_2000_C90
./remap.pl -par regrida2x2.5.par -in jan2010_o3_shindell_144x90x49x12_April1850 -out jan2010_o3_shindell_144x90x49x12_April1850_C90

#-- CH4 Natural sources - have been adjusted and rescaled during pre-processing step in /discover/nobackup/gfaluveg/PRE/make_1x1_nat_sources/CH4

./remap.pl -par regridCH4natural.par -in CH4WETL+TUNDRA_1X1_temp -out CH4WETL+TUNDRA_C90_Dec_2009
./remap.pl -par regridCH4natural.par -in CH4TRMITE_1X1_temp -out CH4TRMITE_C90_Dec_2009
./remap.pl -par regridCH4natural.par -in CH4SOILABS_1X1_temp -out CH4SOILABS_C90_Dec_2009

#-- NOx natural sources - have been adjusted and rescaled during pre-processing step in /discover/nobackup/gfaluveg/PRE/make_1x1_nat_sources/NOx
./remap.pl -par regridNOxnatural.par -in NOx_Soil_GEIA_1x1_half -out NOx_Soil_GEIA_C90_Dec_2009

#-- convert to GISS and add headers
./add_header2.ksh NOx_soil_AR5_1850_C90 m NOx soil C90 N 1850
./add_header2.ksh CH4_termites_AR5_1850_C90 a CH4 termites C90 N 1850
./add_header2.ksh CH4_soil_absorption_AR5_1850_C90 a CH4 soil_absorption C90 N 1850
./add_header2.ksh CH4_Wetlands_and_Tundra_AR5_1850_C90 m CH4 Wetlands_and_Tundra C90 N 1850

#-- Isoprene and terpenes from vegetation - titles must be modified using change_title.f (not performed here)
#./remap.pl -par regridNat.par -in ORCHIDEE_Isoprene_1990_4x5_h_2x2.5gf_h -out ORCHIDEE_Isoprene_1990_C90from2x2.5_h
#./remap.pl -par regridNat.par -in ORCHIDEE_Terpenes_1990_4x5_h_2x2.5gf_h -out ORCHIDEE_Terpenes_1990_C90from2x2.5_h
#./remap.pl -par regridNat.par -in ORCHIDEE_ORVOC_1990_4x5_h_2x2.5gf_h -out ORCHIDEE_ORVOC_1990_C90from2x2.5_h
./remap.pl -par ncregrid-up.par -in bvoc_pr.nc -out bvoc_pr_C90.nc

#-- Alkenes and Paraffin from vegetation (*forC90 are obtained using /gpfsm/dnb53/gfaluveg/PRE/make_F_nat_sources/VOCs/convert1x1.f90)
./remap.pl -par regridb1x1.par -in Alkenes_vegetation_GEIA_1x1_forC90 -out Alkenes_vegetation_GEIA_1x1_C90_h
./remap.pl -par regridb1x1.par -in Paraffin_vegetation_GEIA_1x1_forC90 -out Paraffin_vegetation_GEIA_1x1_C90_h

#-- volcanoes
./remap.pl -par ncregrid-ij.par -in SO2_volc_2000_1x1_AEROCOM.nc -out SO2_volc_2000_AEROCOM_C90.nc

#-- NH3 oceanic source
#./remap.pl -par regridNH3oc.par -in NH3hCON_OCEAN_Apr09_1x1_h -out NH3hCON_OCEAN_Apr09_C90_h
./remap.pl -par regridNH3oc.par -in NH3hCON_OCEANflux_Jan10_1x1_h -out NH3hCON_OCEANflux_Jan10_C90_h 

#-- DMS water conc. Kettle & Andreae 1996
./remap.pl -par regrida1x1.par -in DMS_Kettle_Andeae_1x1 -out DMS_Kettle_Andeae_C90

#-- Terpene from Guenther
./remap.pl -par regrida1x1.par -in terp_Guenther_1x1 -out terp_Guenther_C90

#-- SULFATE_SA, DMS_FIELD, SO2_FIELD
./remap.pl -par regridM.par -in sulfate_fakeM23_M_SA -out sulfate_fakeM23_M_SA_C90
./remap.pl -par regridM.par -in sulfate_pi_fakeM23_M_SA -out sulfate_pi_fakeM23_M_SA_C90
./remap.pl -par regrida2x2.5.par -in sulfate_fakeM23_M_SA_2x2.5gf -out sulfate_fakeM23_M_SA_2x2.5gf_C90
./remap.pl -par regrida2x2.5.par -in sulfate_pi_fakeM23_M_SA_2x2.5gf -out sulfate_pi_fakeM23_M_SA_2x2.5gf_C90
./remap.pl -par regridM.par -in dms_conc -out dms_conc_C90
./remap.pl -par regridM.par -in so2_conc -out so2_conc_C90
./remap.pl -par regrida2x2.5.par -in so2_conc_2x2.5gf -out so2_conc_2x2.5gf_C90
./remap.pl -par regrida2x2.5.par -in dms_conc_2x2.5gf -out dms_conc_2x2.5gf_C90

#-- Dust tracers
#ERS
#./remap.pl -par regridd2x2.5.par -in ERS1_1993_MONTHLY.72x46.threshold-13 -out ERS1_1993_MONTHLY.72x46.threshold-13_C90
./remap.pl -par regrid0.5x0.5.par -in ERS1_1993_MONTHLY.720x360.threshold-13 -out ERS1_1993_MONTHLY.720x360.threshold-13_C90

#GIN
#./remap.pl -par regridb2x2.5.par -in Ginoux2001_source_VegMask_144x90 -out Ginoux2001_source_VegMask_144x90_C90
#./remap.pl -par regridb2x2.5.par -in Ginoux_source_v2009_VegMask_144x90 -out Ginoux_source_v2009_VegMask_144x90_C90
#./remap.pl -par regridb2x2.5.par -in Ginoux_source_v2009_NoVegMask_144x90 -out Ginoux_source_v2009_NoVegMask_144x90_C90
./remap.pl -par regrid0.5x0.5.par -in Ginoux_2001_updated_hires_withVegMask_v2009_0.5x0.5 -out Ginoux_2001_updated_hires_withVegMask_v2009_C90
#-- CO
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_CO_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_CO_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_CO_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_CO_ships_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_CO_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_CO_decadalmonthlymean1850_C90_Dec_2009.nc

#-- CH4
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_CH4_anthropogenic_1850_0.5x0.5_v0_14_10_2009.nc -out IPCC_emissions_CH4_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_CH4_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_CH4_ships_1850_C90_Dec_2009.nc  
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_CH4_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_CH4_decadalmonthlymean1850_C90_Dec_2009.nc 

#-- NH3
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_NH3_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_NH3_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_NH3_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_NH3_decadalmonthlymean1850_C90_Dec_2009.nc 

#-- BCII
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_BCII_anthropogenic_1850_0.5x0.5_v1_21_05_2009.nc -out IPCC_emissions_BCII_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_BCII_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_BCII_ships_1850_C90_Dec_2009.nc  

#-- OCII
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_OCII_anthropogenic_1850_0.5x0.5_v1_21_05_2009.nc -out IPCC_emissions_OCII_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_OCII_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_OCII_ships_1850_C90_Dec_2009.nc  

#-- SO2
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_SO2_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_SO2_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_SO2_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_SO2_ships_1850_C90_Dec_2009.nc  
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_SO2_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_SO2_decadalmonthlymean1850_C90_Dec_2009.nc 

#-- BCB
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_BCB_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_BCB_decadalmonthlymean1850_C90_Dec_2009.nc 

#-- OCB
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_OCB_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_OCB_decadalmonthlymean1850_C90_Dec_2009.nc 



#-- NOx
# first perform NOx=NO*14./30.
ncflint -c -v emiss_agr,emiss_awb,emiss_dom,emiss_ene,emiss_ind,emiss_tra,emiss_wst -w .46666666666666666666,0.0 IPCC_emissions_NO_anthropogenic_1850_0.5x0.5_v1_07_05_2009.nc IPCC_emissions_NO_anthropogenic_1850_0.5x0.5_v1_07_05_2009.nc -o NOx_anthropo.nc
ncflint -c -v emiss_shp -w .46666666666666666666,0.0 IPCC_emissions_NO_ships_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_NO_ships_1850_0.5x0.5_v1_20_04_2009.nc -o NOx_ships.nc
ncflint -c -v grassfire,forestfire -w .46666666666666666666,0.0  IPCC_GriddedBiomassBurningEmissions_NOx_decadalmonthlymean1850_v1.nc IPCC_GriddedBiomassBurningEmissions_NOx_decadalmonthlymean1850_v1.nc -o NOx_Biomass.nc
./remap.pl -par ncregrid-ijl.par -in NOx_anthropo.nc -out IPCC_emissions_NOx_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in NOx_ships.nc -out IPCC_emissions_NOx_ships_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in NOx_Biomass.nc -out IPCC_GriddedBiomassBurningEmissions_NOx_decadalmonthlymean1850_C90_Dec_2009.nc

#-- Alkenes
# first perform Alkenes = propene/42.0 + other_alkenes_and_alkynes/67.0 + ethene/28.0
ncflint -c -v emiss_shp -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_ships_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_ships_1850_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_ships.nc 
ncflint -c -v emiss_shp -w 1.0,0.03571428571428571428 propene+oalkenes_ships.nc IPCC_emissions_ethene_ships_1850_0.5x0.5_v1_20_04_2009.nc -o alkenes_ships.nc 

ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_anthropogenic.nc
ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 1.0,0.03571428571428571428 propene+oalkenes_anthropogenic.nc IPCC_emissions_ethene_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -o alkenes_anthropogenic.nc

ncflint -c -v grassfire,forestfire -w 0.02380952380952380952,0.01492537313432835820 IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean1850_v1.nc IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean1850_v1.nc -o propene+oalkenes_Biomass.nc 
ncflint -c -v grassfire,forestfire -w 1.0,0.03571428571428571428 propene+oalkenes_Biomass.nc IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean1850_v1.nc -o alkenes_Biomass.nc

./remap.pl -par ncregrid-ijl.par -in alkenes_anthropogenic.nc -out IPCC_emissions_alkenes_anthropogenic_1850_C90_Dec_2009.nc

./remap.pl -par ncregrid-ijl.par -in alkenes_ships.nc -out IPCC_emissions_alkenes_ships_1850_C90_Dec_2009.nc

./remap.pl -par ncregrid-ijl.par -in alkenes_Biomass.nc -out IPCC_GriddedBiomassBurningEmissions_alkenes_decadalmonthlymean1850_C90_Dec_2009.nc

#-- Paraffin
# first perform Paraffin = propane/44.0 + pentanes/72.0 + butanes/57.8 + hexanes_and_higher_alkanes/106.8 + ethane/30.0 + ketones/75.3
ncflint -c -v emiss_shp -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_ships_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_ships_1850_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_ships.nc
ncflint -c -v emiss_shp -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_ships_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_ships_1850_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_ships.nc
ncflint -c -v emiss_shp -w 0.03333333333333333333,0.0 IPCC_emissions_ethane_ships_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ethane_ships_1850_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_ships.nc  #- note the weight == 0.0 there is no ketones_ships file
ncflint -c -v emiss_shp -w 1.0,1.0 propane+pentanes_ships.nc butanes+hexanes_ships.nc -o propane+pentanes+butanes+hexanes_ships.nc
ncflint -c -v emiss_shp -w 1.0,1.0 propane+pentanes+butanes+hexanes_ships.nc ethane+ketones_ships.nc -o paraffin_ships.nc

ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_anthropogenic.nc
ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_anthropogenic.nc
ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 0.03333333333333333333,0.01328021248339973439 IPCC_emissions_ethane_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ketones_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_anthropogenic.nc
ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 1.0,1.0 propane+pentanes_anthropogenic.nc butanes+hexanes_anthropogenic.nc -o propane+pentanes+butanes+hexanes_anthropogenic.nc
ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 1.0,1.0 propane+pentanes+butanes+hexanes_anthropogenic.nc ethane+ketones_anthropogenic.nc -o paraffin_anthropogenic.nc

ncflint -c -v grassfire,forestfire -w 0.02272727272727272727,0.01388888888888888888 IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean1850_v1.nc IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean1850_v1.nc -o propane+pentanes_Biomass.nc
ncflint -c -v grassfire,forestfire -w 0.01730103806228373702,0.00936329588014981273 IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean1850_v1.nc IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean1850_v1.nc -o butanes+hexanes_Biomass.nc
ncflint -c -v grassfire,forestfire -w 0.03333333333333333333,0.01328021248339973439 IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean1850_v1.nc IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean1850_v1.nc -o ethane+ketones_Biomass.nc
ncflint -c -v grassfire,forestfire -w 1.0,1.0 propane+pentanes_Biomass.nc butanes+hexanes_Biomass.nc -o propane+pentanes+butanes+hexanes_Biomass.nc
ncflint -c -v grassfire,forestfire -w 1.0,1.0 propane+pentanes+butanes+hexanes_Biomass.nc ethane+ketones_Biomass.nc -o paraffin_Biomass.nc

# then do regridding
./remap.pl -par ncregrid-ijl.par -in paraffin_ships.nc -out IPCC_emissions_paraffin_ships_1850_C90_Dec_2009.nc

./remap.pl -par ncregrid-ijl.par -in paraffin_anthropogenic.nc -out IPCC_emissions_paraffin_anthropogenic_1850_C90_Dec_2009.nc

./remap.pl -par ncregrid-ijl.par -in paraffin_Biomass.nc -out IPCC_GriddedBiomassBurningEmissions_paraffin_decadalmonthlymean1850_C90_Dec_2009.nc

#-- convert to GISS format and add headers
    module avail tool/idl
    module load tool/idl-6.4
    ulimit -s 6000000
    ulimit -v unlimited
    idl ./nc2giss.pro    # terpene isoprene
    cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro
    ./run_1850_anthro_CS.ksh
    ./run_1850_ships_CS.ksh
    ./run_1850_BiomassBurning_CS.ksh
    ./run_1850_anthro_alkenes_paraffin_CS.ksh
#    ./zero_veg_CS.ksh
    ./zero_aircraft_CS.ksh 

    cp HTAP_surface_co_2001_aggregated_C90.nc /archive/u/dgueyffi/codirect_for_ACCMIP/
    cd /archive/u/dgueyffi/codirect_for_ACCMIP
    idl convert_co_htap_for_accmip_C90.pro 

#-- transient files
cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro/out_auto_C90
./make_surface_transient_Alkenes_C90.sh
./make_surface_transient_Paraffin_C90.sh
./make_surface_transient_NOx_C90.sh
./make_surface_transient_CO_C90.sh
./make_surface_transient_CH4_C90.sh
./make_surface_transient_BCB.sh
./make_surface_transient_BCII.sh
./make_surface_transient_NH3.sh
./make_surface_transient_OCB.sh
./make_surface_transient_OCII.sh
./make_surface_transient_SO2.sh
./loop_headers.ksh


#-- check that global integrals  are conserved through regridding
    cd /gpfsm/dnb53/gfaluveg/AR5_emissions/v5_anthro/
idl -e "testint1r8" -args gsin/N2O_IC_M23_4x5_6.17_conc_2x2.5_conc AR5_C90/N2O_IC_from_2x2.5_C90_Dec_2009
idl -e "testint1r8" -args gsin/CFC_IC_M23_4x5_6.17_conc_2x2.5_conc AR5_C90/CFC_IC_from_2x2.5_C90_Dec_2009
idl -e "testint1r8" -args gsin/CH4_IC_M23_4x5_6.17_conc_2x2.5_conc AR5_C90/CH4_IC_from_2x2.5_C90_Dec_2009
idl -e "testint1r8" -args gsin/Ox_init_cond_M23_4x5_conc_2x2.5_conc AR5_C90/Ox_IC_from_2x2.5_C90_Dec_2009
idl -e "testint1r8" -args gsin/CO_init_cond_M23_conc_2x2.5_conc AR5_C90/CO_IC_from_2x2.5_C90_Dec_2009
idl -e "testint" -args temp_2x2.5/sulfate_pi_fakeM23_M_SA_2x2.5gf AR5_C90/sulfate_pi_fakeM23_M_SA_2x2.5gf_C90
idl -e "testint" -args temp_2x2.5/dms_conc_2x2.5gf AR5_C90/dms_conc_2x2.5gf_C90
idl -e "testint" -args temp_2x2.5/so2_conc_2x2.5gf AR5_C90/so2_conc_2x2.5gf_C90 
idl -e "testint0" -args ERS1_1993_MONTHLY.144x90.threshold-13 AR5_C90/ERS1_1993_MONTHLY.720x360.threshold-13_C90
idl -e "testint0" -args Ginoux_source_v2009_VegMask_144x90 AR5_C90/Ginoux_2001_updated_hires_withVegMask_v2009_C90
idl -e "testint" -args AR5_emis/F/T/CO_ind_AR5_1850-2000_2x2.5_h AR5_C90/transient/CO_ind_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/CO_tra_AR5_1850-2000_2x2.5_h AR5_C90/transient/CO_tra_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/CO_wst_AR5_1850-2000_2x2.5_h AR5_C90/transient/CO_wst_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/CO_awb_AR5_1850-2000_2x2.5_h AR5_C90/transient/CO_awb_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/CO_dom_AR5_1850-2000_2x2.5_h AR5_C90/transient/CO_dom_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/CO_forestfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/CO_forestfire_AR5_1900-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/CO_grassfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/CO_grassfire_AR5_1900-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/m_CO_shp_AR5_1850-2000_2x2.5_h AR5_C90/transient/CO_shp_AR5_1850-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/CO_slv_AR5_1990-2000_2x2.5_h AR5_C90/transient/CO_slv_AR5_1990-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/CO_ene_AR5_1850-2000_2x2.5_h AR5_C90/transient/CO_ene_AR5_1850-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/CO_agr_AR5_1990-2000_2x2.5_h AR5_C90/transient/CO_agr_AR5_1990-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/NOx_air_AR5_1910-2000_2x2.5 AR5_C90/NOx_air_AR5_1910-2000_C90
idl -e "testint" -args AR5_emis/F/NAT/NOx_Soil_GEIA_2x2.5_HALF_h AR5_C90/NOx_soil_AR5_1850_C90_h
idl -e "testint" -args AR5_emis/F/T/NOx_awb_AR5_1850-2000_2x2.5_h AR5_C90/transient/NOx_awb_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/NOx_dom_AR5_1850-2000_2x2.5_h AR5_C90/transient/NOx_dom_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/NOx_ene_AR5_1850-2000_2x2.5_h AR5_C90/transient/NOx_ene_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/NOx_forestfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/NOx_forestfire_AR5_1900-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/NOx_grassfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/NOx_grassfire_AR5_1900-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/NOx_ind_AR5_1850-2000_2x2.5_h AR5_C90/transient/NOx_ind_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/m_NOx_shp_AR5_1850-2000_2x2.5_h AR5_C90/transient/NOx_shp_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/NOx_tra_AR5_1850-2000_2x2.5_h AR5_C90/transient/NOx_tra_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/NOx_wst_AR5_1850-2000_2x2.5_h AR5_C90/transient/NOx_wst_AR5_1850-2000_C90_h
idl -e "testint4" -args AR5_emis/F/T/NOx_agr_AR5_1850-2000_2x2.5_h AR5_C90/transient/NOx_agr_AR5_1850-2000_C90_h
#idl -e "testint" -args ORCHIDEE_Terpenes_1990_2x2.5_h 
#idl -e "testint" -args ORCHIDEE_ORVOC_1990_2x2.5_h 
idl -e "testint" -args AR5_emis/F/NAT/Alkenes_vegetation_GEIA_2x2.5_h_1 AR5_C90/Alkenes_vegetation_GEIA_1x1_C90_h
idl -e "testint" -args AR5_emis/F/T/m_Alkenes_shp_AR5_1850-2000_2x2.5_h AR5_C90/transient/Alkenes_shp_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/Alkenes_wst_AR5_1850-2000_2x2.5_h AR5_C90/transient/Alkenes_wst_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/Alkenes_dom_AR5_1850-2000_2x2.5_h AR5_C90/transient/Alkenes_dom_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/Alkenes_forestfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/Alkenes_forestfire_AR5_1900-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/Alkenes_grassfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/Alkenes_grassfire_AR5_1900-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/Alkenes_ind_AR5_1850-2000_2x2.5_h AR5_C90/transient/Alkenes_ind_AR5_1850-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/Alkenes_tra_AR5_1850-2000_2x2.5_h AR5_C90/transient/Alkenes_tra_AR5_1850-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/Alkenes_ene_AR5_1850-2000_2x2.5_h AR5_C90/transient/Alkenes_ene_AR5_1850-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/Alkenes_awb_AR5_1890-2000_2x2.5_h AR5_C90/transient/Alkenes_awb_AR5_1890-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/Alkenes_agr_AR5_1990-2000_2x2.5_h AR5_C90/transient/Alkenes_agr_AR5_1990-2000_C90_h
idl -e "testint" -args AR5_emis/F/NAT/Paraffin_vegetation_GEIA_2x2.5_h_1 AR5_C90/Paraffin_vegetation_GEIA_1x1_C90_h
idl -e "testint" -args AR5_emis/F/T/m_Paraffin_shp_AR5_1850-2000_2x2.5_h AR5_C90/transient/Paraffin_shp_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/Paraffin_wst_AR5_1850-2000_2x2.5_h AR5_C90/transient/Paraffin_wst_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/Paraffin_dom_AR5_1850-2000_2x2.5_h AR5_C90/transient/Paraffin_dom_AR5_1850-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/Paraffin_forestfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/Paraffin_forestfire_AR5_1900-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/Paraffin_grassfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/Paraffin_grassfire_AR5_1900-2000_C90_h
idl -e "testint" -args AR5_emis/F/T/Paraffin_ind_AR5_1850-2000_2x2.5_h AR5_C90/transient/Paraffin_ind_AR5_1850-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/Paraffin_tra_AR5_1850-2000_2x2.5_h AR5_C90/transient/Paraffin_tra_AR5_1850-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/Paraffin_slv_AR5_1860-2000_2x2.5_h AR5_C90/transient/Paraffin_slv_AR5_1860-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/Paraffin_ene_AR5_1850-2000_2x2.5_h AR5_C90/transient/Paraffin_ene_AR5_1850-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/Paraffin_awb_AR5_1890-2000_2x2.5_h AR5_C90/transient/Paraffin_awb_AR5_1890-2000_C90_h
idl -e "testint3" -args AR5_emis/F/T/Paraffin_agr_AR5_1990-2000_2x2.5_h AR5_C90/transient/Paraffin_agr_AR5_1990-2000_C90_h
idl -e "testint" -args AR5_emis/F/HTAP_codirect_emissions_2x2.5_h AR5_C90/HTAP_codirect_emissions_C90_h
idl -e "testint" -args NH3hCON_OCEANflux_Jan10_2x2.5_h AR5_C90/NH3hCON_OCEANflux_Jan10_C90title_h
idl -e "testint" -args NH3_forestfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/NH3_forestfire_AR5_1900-2000_C90_h
idl -e "testint" -args NH3_grassfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/NH3_grassfire_AR5_1900-2000_C90_h
idl -e "testint" -args NH3_agr_AR5_1850-2000_2x2.5_h AR5_C90/transient/NH3_agr_AR5_1850-2000_C90_h
idl -e "testint" -args NH3_awb_AR5_1850-2000_2x2.5_h AR5_C90/transient/NH3_awb_AR5_1850-2000_C90_h
idl -e "testint" -args NH3_dom_AR5_1850-2000_2x2.5_h AR5_C90/transient/NH3_dom_AR5_1850-2000_C90_h
idl -e "testint" -args NH3_ind_AR5_1850-2000_2x2.5_h AR5_C90/transient/NH3_ind_AR5_1850-2000_C90_h
idl -e "testint4" -args NH3_ene_AR5_1850-2000_2x2.5_h AR5_C90/transient/NH3_ene_AR5_1850-2000_C90_h
idl -e "testint4" -args NH3_tra_AR5_1850-2000_2x2.5_h AR5_C90/transient/NH3_tra_AR5_1850-2000_C90_h
#idl -e "testint" -args SO2_volc_conti2000_HR2x2.5.AEROCOM 
idl -e "testint" -args DMS_Kettle_Andreae_2x2.5 AR5_C90/DMS_Kettle_Andeae_C90
idl -e "testint" -args BCII_awb_AR5_1850-2000_2x2.5_h AR5_C90/transient/BCII_awb_AR5_1850-2000_C90_h
idl -e "testint" -args BCII_dom_AR5_1850-2000_2x2.5_h AR5_C90/transient/BCII_dom_AR5_1850-2000_C90_h
idl -e "testint3" -args BCII_ene_AR5_1850-2000_2x2.5_h AR5_C90/transient/BCII_ene_AR5_1850-2000_C90_h
idl -e "testint" -args BCII_ind_AR5_1850-2000_2x2.5_h AR5_C90/transient/BCII_ind_AR5_1850-2000_C90_h
idl -e "testint" -args BCII_tra_AR5_1850-2000_2x2.5_h AR5_C90/transient/BCII_tra_AR5_1850-2000_C90_h
idl -e "testint" -args BCII_wst_AR5_1850-2000_2x2.5_h AR5_C90/transient/BCII_wst_AR5_1850-2000_C90_h
idl -e "testint" -args BCII_shp_AR5_1850-2000_2x2.5_h AR5_C90/transient/BCII_shp_AR5_1850-2000_C90_h
idl -e "testint" -args OCII_awb_AR5_1850-2000_2x2.5_h AR5_C90/transient/OCII_awb_AR5_1850-2000_C90_h
idl -e "testint" -args OCII_dom_AR5_1850-2000_2x2.5_h AR5_C90/transient/OCII_dom_AR5_1850-2000_C90_h
idl -e "testint3" -args OCII_ene_AR5_1850-2000_2x2.5_h AR5_C90/transient/OCII_ene_AR5_1850-2000_C90_h
idl -e "testint" -args OCII_ind_AR5_1850-2000_2x2.5_h AR5_C90/transient/OCII_ind_AR5_1850-2000_C90_h
idl -e "testint" -args OCII_tra_AR5_1850-2000_2x2.5_h AR5_C90/transient/OCII_tra_AR5_1850-2000_C90_h
idl -e "testint" -args OCII_wst_AR5_1850-2000_2x2.5_h AR5_C90/transient/OCII_wst_AR5_1850-2000_C90_h
idl -e "testint" -args OCII_shp_AR5_1850-2000_2x2.5_h AR5_C90/transient/OCII_shp_AR5_1850-2000_C90_h
idl -e "testint" -args SO2_shp_AR5_1850-2000_2x2.5_h AR5_C90/transient/SO2_shp_AR5_1850-2000_C90_h
idl -e "testint" -args SO2_awb_AR5_1850-2000_2x2.5_h AR5_C90/transient/SO2_awb_AR5_1850-2000_C90_h
idl -e "testint" -args SO2_dom_AR5_1850-2000_2x2.5_h AR5_C90/transient/SO2_dom_AR5_1850-2000_C90_h
idl -e "testint" -args SO2_ind_AR5_1850-2000_2x2.5_h AR5_C90/transient/SO2_ind_AR5_1850-2000_C90_h
idl -e "testint" -args SO2_tra_AR5_1850-2000_2x2.5_h AR5_C90/transient/SO2_tra_AR5_1850-2000_C90_h
idl -e "testint" -args SO2_wst_AR5_1850-2000_2x2.5_h AR5_C90/transient/SO2_wst_AR5_1850-2000_C90_h
idl -e "testint3" -args SO2_ene_AR5_1850-2000_2x2.5_h AR5_C90/transient/SO2_ene_AR5_1850-2000_C90_h
idl -e "testint" -args BCB_forestfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/BCB_forestfire_AR5_1900-2000_C90_h
idl -e "testint" -args BCB_grassfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/BCB_grassfire_AR5_1900-2000_C90_h
idl -e "testint" -args OCB_forestfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/OCB_forestfire_AR5_1900-2000_C90_h
idl -e "testint" -args OCB_grassfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/OCB_grassfire_AR5_1900-2000_C90_h
idl -e "testint" -args SO2_forestfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/SO2_forestfire_AR5_1900-2000_C90_h
idl -e "testint" -args SO2_grassfire_AR5_1900-2000_2x2.5_h AR5_C90/transient/SO2_grassfire_AR5_1900-2000_C90_h
