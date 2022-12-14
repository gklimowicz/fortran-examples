#!/usr/bin/perl

# creates one 2D map  at specified level

###### ---------Script--------- ########
do 'user_input.s';

chdir $DataDir;

print "obs: $ObsDir$ObsFilename \n";
print "variable_obs $variable_obs \n";

@data_array = ("DJF","JJA");
for my $months (@data_array){

  $OutputFileName = "$variable_obs.$months.map_lev$ilev.$RUN2.nc";
  print "$OutputFileName.\n";

  if (defined($depth)){
    print "Extracting depth at level $ilev \n";
    system "ncks -O -F -d $depth,$ilev,$ilev,1 -v $variable_obs $ObsDir$ObsFilename dummy.nc";
  }else{
    system "ncks -O -v $variable_obs $ObsDir$ObsFilename dummy.nc";
  }

  if ($months eq "DJF"){
    system "ncks -O -F -d mon,1,2 -v $variable_obs dummy.nc dummy1.nc";
    system "ncks -A -F -d mon,11,11 -v $variable_obs dummy.nc dummy1.nc";
    system "ncwa -O -v $variable_obs -a mon dummy1.nc $OutputFileName";
  }
  if ($months eq "JJA"){
    system "ncks -O -F -d mon,6,8 -v $variable_obs dummy.nc dummy1.nc";
    system "ncwa -O -a mon -v $variable_obs dummy1.nc $OutputFileName";
}
system "ncatted -O -a _FillValue,$variable_obs,d,f, $OutputFileName";
system "ncatted -O -a _FillValue,$variable_obs,c,f,-1e30 $OutputFileName";   
system "ncatted -O -a units,$lat_obs,c,c,'degrees_north' $OutputFileName";
system "ncatted -O -a units,$lon_obs,c,c,'degrees_east' $OutputFileName";

system "rm -R -f dummy*.nc";

#rename lat_obs/lon_obs to "lat" and "lon"
print "outputfile $OutputFileName";
if ($lat_obs ne "lat") {
    system "ncrename -O -v $lat_obs,lat -d $lat_obs,lat $OutputFileName";
    system "ncrename -O -v $lon_obs,lon -d $lon_obs,lon $OutputFileName";
}

#rename variable_obs to variable
if ($variable_obs ne $variable) {
system "ncrename -O -v $variable_obs,$variable $OutputFileName";
}
}

#####
$months = "ANN";
print "$months \n";

if (index($variable_obs, "_mon") != -1) {
   print "'$variable_obs' contains mon.\n";
  $variable_obs =~ s/mon/ann/g; 
}

$OutputFileName = "$variable_obs.$months.map_lev$ilev.$RUN2.nc";
print "$OutputFileName.\n";

if (defined($depth)){
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -O -F -d $depth,$ilev,$ilev,1 -v $variable_obs $ObsDir$ObsFilename $OutputFileName";
}elsif (defined($mon)){
  system "ncks -O -v $variable_obs $ObsDir$ObsFilename dummy.nc";
  system "ncwa -O -v $variable_obs -a mon dummy.nc dummy1.nc";
  system "ncks -O -x -v mon dummy1.nc $OutputFileName";
}else{
  system "ncks -O -v $variable_obs $ObsDir$ObsFilename $OutputFileName";
}

system "ncatted -O -a _FillValue,$variable_obs,d,f, $OutputFileName";
system "ncatted -O -a _FillValue,$variable_obs,c,f,-1e30 $OutputFileName";   
system "ncatted -O -a units,$lat_obs,c,c,'degrees_north' $OutputFileName";
system "ncatted -O -a units,$lon_obs,c,c,'degrees_east' $OutputFileName";
system "rm -R -f dummy*.nc";


#rename lat_obs/lon_obs to "lat" and "lon"
if ($lat_obs ne "lat") {
    system "ncrename -O -v $lat_obs,lat -d $lat_obs,lat $OutputFileName";
    system "ncrename -O -v $lon_obs,lon -d $lon_obs,lon $OutputFileName";
}

#rename variable_obs to variable
if ($variable_obs ne $variable) {
system "ncrename -O -v $variable_obs,$variable $OutputFileName";
}
