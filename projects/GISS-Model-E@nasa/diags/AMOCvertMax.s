#!/usr/bin/perl

##### -------- Inputs ------- #####

$nctag = ojl;
$variable = sf_Atl;

$lat = lato2;
$latVal = 26.5;
###### ---------Script--------- ########

chdir $DataDir; 

$OutputFileName = "$variable.$yrini-$yrend.AMOCvertMax_lat$latVal.$RUN.nc";
print "Output File $OutputFileName \n";

$year = $yrini;
while ($year <= $yrend)
{
  system "scaleacc ANN$year.acc$RUN.nc $nctag";
  system "ncks -O -a -v $variable ANN$year.$nctag$RUN.nc dummy.nc";
  system "ncks -a -v $variable -d $lat,$latVal dummy.nc dummy1.nc";
  system "ncwa -y max -v $variable dummy1.nc dum_$year.nc";
  system "rm -f dummy*.nc";
  $year=$year + 1;
}

system "ncecat dum_*.nc $OutputFileName";

system "rm -f dum*.nc";

