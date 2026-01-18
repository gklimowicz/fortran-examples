#!/usr/bin/perl

# lat/dep plot from 3D model data form ojl
# in different ocean basins, Pacific and Atlantic

do 'user_input.s';

$variable_moc= "sf";
$nctag = ojl;
$underscore = "_";
$lon = lono2;

###### ---------Script--------- ########

chdir $DataDir; 

$OutputFileName = "$variable_moc.ANN$yrini-$yrend.MOCBasin.$RUN.nc";
print "$OutputFileName\n";

@data_array = ("Atl","Pac","Ind","Glo");
for my $Basin (@data_array){

$variable = "$variable_moc$underscore$Basin";
system "scaleacc ANN$yrini-$yrend.acc$RUN.nc $nctag";
system "ncks -v $variable ANN$yrini-$yrend.$nctag$RUN.nc dummy.nc";
system "ncks -A -a -v $variable dummy.nc $OutputFileName";
system "ncatted -O -a _FillValue,$variable,c,f,-6.849315e+26 $OutputFileName";

system "rm -f dummy*.nc";
}
