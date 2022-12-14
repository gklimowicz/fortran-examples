#!/usr/bin/perl

# creates a lat/dep plot from gridded observations
# in different ocean basins, Pacific and Atlantic

###### ---------Script--------- ########
do 'user_input.s';

chdir $DataDir;

if (index($variable_obs, "_mon") != -1) {
   print "'$variable_obs' contains mon.\n";
  $variable_obs =~ s/mon/ann/g; 
}

$ann="ann";
if ($variable_obs eq "flux_molCm2yr_tak") {
    $variable_obs = "$variable_obs$underscore$ann";
}  
 
print "variable_obs $variable_obs \n";
print "obsfilename $ObsFilename\n";
print "extract the variable_obs \n";
system "ncks -O -v $variable_obs $ObsDir$ObsFilename dummy.nc";

$OutputFileName = "$variable_obs.Basin.$RUN2.nc";
print "outputfilename $OutputFileName\n";


if ($nctag eq "taij"){
 print "variable is flux \n";
 $lat_obs = lata;
 $lon_obs = lona;
}


###### ---------Atlantic Basin

$basin = "Atl";
print "$basin\n";

system "ncrename -v $variable_obs,AtlMask dummy.nc";

if ($RUN2 eq "Takahashi_onEgrid"){
  print "Rename the lat/lon in dummy.nc\n";
  system "ncrename -O -v $lat_obs,lat -d $lat_obs,lat dummy.nc";
  system "ncrename -O -v $lon_obs,lon -d $lon_obs,lon dummy.nc";
}

if ($area eq "oxyp"){
  $mask = "AtlMaskO.nc";
  system "ncbo --op_typ=multiply $ObsDir$mask dummy.nc dummy1.nc";
  if (! defined($depth)){
  system "ncks -O -F -d zoc,1,1,1 -v AtlMask dummy1.nc dummy1.nc";
  }
}else{
  $mask = "AtlMaskA.nc";
  system "ncbo --op_typ=multiply $ObsDir$mask dummy.nc dummy1.nc";
}

system "ncwa -v AtlMask -a lon dummy1.nc dummy2.nc";
$variable_obs_new = "$variable_obs$basin";
print "$variable_obs_new \n";
system "ncrename -v AtlMask,$variable_obs_new dummy2.nc dummy3.nc";
# Remove unwanted terms
system "ncks -O -x -v lon dummy3.nc $OutputFileName";

# flip depth DOESN't Work...
#system "ncpdq -O -h -a -zoc dummy3.nc $OutputFileName";
system "ncatted -O -a _FillValue,$variable_obs_new,d,f, $OutputFileName"; 
system "ncatted -O -a _FillValue,$variable_obs_new,c,f,-1e30 $OutputFileName"; 
system "rm -R -f dummy*";

###### ---------Pacific Basin
$basin = "Pac";
print "$basin\n";

print "extract the variable_obs \n";
system "ncks -O -v $variable_obs $ObsDir$ObsFilename dummy.nc";
print "$ObsFilename\n";
system "ncrename -v $variable_obs,PacMask dummy.nc";

if ($RUN2 eq "Takahashi_onEgrid"){
  print "Rename the lat/lon in dummy.nc\n";
  system "ncrename -O -v $lat_obs,lat -d $lat_obs,lat dummy.nc";
  system "ncrename -O -v $lon_obs,lon -d $lon_obs,lon dummy.nc";
}

if ($area eq "oxyp"){
  $mask = "PacMaskO.nc";
  system "ncbo --op_typ=multiply $ObsDir$mask dummy.nc dummy1.nc";
  if (! defined($depth)){
  system "ncks -O -F -d zoc,1,1,1 -v PacMask dummy1.nc dummy1.nc";
}
}else{
  $mask = "PacMaskA.nc";
  system "ncbo --op_typ=multiply $ObsDir$mask dummy.nc dummy1.nc";
}

system "ncwa -v PacMask -a lon dummy1.nc dummy2.nc";
$ann="ann";
$variable_obs_new = "$variable_obs$basin";
print "$variable_obs_new \n";
system "ncrename -v PacMask,$variable_obs_new dummy2.nc";
# Remove unwanted terms
system "ncks -O -x -v lon dummy2.nc dummyPac.nc";

system "ncks -A -v $variable_obs_new dummyPac.nc $OutputFileName";
system "ncatted -O -a _FillValue,$variable_obs_new,d,f, $OutputFileName"; 
system "ncatted -O -a _FillValue,$variable_obs_new,c,f,-1e30 $OutputFileName"; 
system "rm -R -f dummy*";
