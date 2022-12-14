#!/usr/bin/perl

# computess the global average seasonal cycle

###### ---------Script--------- ########
print "Doing SeaonsalCycleObs. Calculates the lat/lon average at a specifc depth per month. \n";
do 'user_input.s';

if ($nctag eq "taij"){
print ("We use flux in mol, (C/m2/yr)\n");
}

$computes = "seasonalCycle_ts";
$new_variablename = "$variable$underscore$computes";
chdir $DataDir;
$OutputFileName = "$variable_obs.$computes.lev$ilev.$ObsFilename";
#substr($OutputFileName,rindex $OutputFileName,'.nc') = '';
#print "HERE $OutputFileName \n";

system "ncks -O -v $area $area$underscore$RUN.nc dummy1.nc";
if ($lat ne $lat_obs){
system "ncrename -O -h -d $lat,$lat_obs -v $lat,$lat_obs dummy1.nc";
system "ncrename -O -h -d $lon,$lon_obs -v $lon,$lon_obs dummy1.nc";
}

if (defined($depth)){
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -A -F -d $depth,$ilev,$ilev,1 -v $variable_obs $ObsDir$ObsFilename dummy1.nc";
}else{
  system "ncks -A -v $variable_obs $ObsDir$ObsFilename dummy1.nc";
}

# create record dimension to match model
system "ncecat -O -h dummy1.nc dummy2.nc";
system "ncpdq -O -h -a mon,record dummy2.nc dummy2.nc";
system "ncwa -O -h -a record dummy2.nc dummy2.nc";

system "ncwa -h -O -v $variable_obs -w $area -a $lat_obs,$lon_obs dummy2.nc $OutputFileName";

# Remove unwanted terms
system "ncks -O -x -v $lon_obs $OutputFileName $OutputFileName";
system "ncks -O -x -v $lat_obs $OutputFileName $OutputFileName";

# Average over depth so we can remove it
if (defined($depth)){
system "ncwa -O -v $variable_obs -a zoc $OutputFileName $OutputFileName";
system "ncks -O -x -v zoc $OutputFileName $OutputFileName";
}

print("$variable_obs \n");
print("$new_variablename \n");
print("$OutputFileName \n");
system "ncrename -O -h -v $variable_obs,$new_variablename $OutputFileName";
system "ncatted -O -a _FillValue,$new_variablename,d,f, $OutputFileName";

# Computes the northern and southern hermisphere averages
if ($variable eq "oicefr"){
   $NH="NortH";
   $SH="SoutH";
   $OutputFileName_noext = "$variable_obs.$computes.lev$ilev.$ObsFilename";
   substr($OutputFileName_noext,rindex $OutputFileName_noext,'.nc') = '';
   
   system "ncwa -h -O -w $area -a $lat,$lon --mask_condition '$lat>0' dummy1.nc dummy3.nc";
   system "ncrename -v $variable,$new_variablename$underscore$NH dummy3.nc";
   $OutputFileName ="$OutputFileName_noext.$NH.nc";
   print "$OutputFileName \n";
   system "ncks -h dummy3.nc $OutputFileName";
   system "ncwa -h -O -w $area -a $lat,$lon --mask_condition '$lat<0' dummy1.nc dummy4.nc";
   $OutputFileName ="$OutputFileName_noext.$SH.nc";
   system "ncrename -v $variable,$new_variablename$underscore$SH dummy4.nc";
   system "ncks -h dummy4.nc $OutputFileName";
   print "$OutputFileName \n";
}

system "rm -f dummy*.nc";


