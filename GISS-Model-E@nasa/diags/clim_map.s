#!/usr/bin/perl

# creates one 2D map climatology at specified level
# for annual/DJF/JJA data

###### ---------Script--------- ########
do 'user_input.s';

print "Doing clim_maps to make a climatology 2D map (lat/lon) \n";
chdir $DataDir;
@data_array = ("ANN","DJF","JJA");

#NOTE: DJF starts at yrini+1;

for my $months (@data_array){
$yrin = $yrini;
if ($months eq "DJF"){
  $yrin = $yrini+1;
}

$InputFileName = "$months$yrin-$yrend.acc$RUN.nc";
print "$InputFileName \n";
 
system "scaleacc $DataDir$InputFileName $nctag";

$OutputFileName = "$variable.$months$yrin-$yrend.map_lev$ilev.$RUN.nc";
print "$OutputFileName.\n";

if ($months eq "ANN"){
  $OutputFileName1=$OutputFileName;
}

if (defined($depth)){
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -O -F -d $depth,$ilev,$ilev,1 -v $variable $months$yrin-$yrend.$nctag$RUN.nc $OutputFileName";
}else{
  system "ncks -O -v $variable $months$yrin-$yrend.$nctag$RUN.nc $OutputFileName";
}

#rename lat/lon from lato/lono
if ($lat ne "lat") {
    system "ncrename -O -v $lat,lat -d $lat,lat $OutputFileName";
    system "ncrename -O -v $lon,lon -d $lon,lon $OutputFileName";
}
}

if ($variable eq "oij_pCO2"){
  system "ncrename -v $variable,LandMask $OutputFileName dummy.nc";
  $mask = "LandMask2D.nc";
  system "ncbo --op_typ='+' $ObsDir$mask dummy.nc dummy1.nc";
  system "ncrename -O -v LandMask,$variable dummy1.nc $OutputFileName";
  system "ncatted -O -a _FillValue,$variable,c,f,-1e30 $OutputFileName";   
  system "ncatted -O -a units,$lat,c,c,'degrees_north' $OutputFileName";
  system "ncatted -O -a units,$lon,c,c,'degrees_east' $OutputFileName";
  system "rm -R -f dummy*";
}

