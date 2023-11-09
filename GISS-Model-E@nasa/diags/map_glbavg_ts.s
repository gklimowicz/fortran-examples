#!/usr/bin/perl

# Computes the global average timeseries at a specified level for 3D field
# by summing the lat/lon and extracting the specified level

###### ---------Script--------- ########

do 'user_input.s';

chdir $DataDir;
$OutputFileName = "$variable.$yrini-$yrend.glbAvg_lev$ilev.$RUN.nc";
print "Output Name $OutputFileName\n";

$year = $yrini;
while ($year <= $yrend) 
{
  system "scaleacc ANN$year.acc$RUN.nc $nctag";
  system "ncks -O -v $variable ANN$year.$nctag$RUN.nc dummy_$year.nc";
  $year=$year + 1;
}

system "ncecat dummy_*.nc $variable.$yrini-$yrend.$RUN.nc";
system "rm -f dummy_*.nc";

#Rename missing_value to FillValue before computation
system "ncatted -O -a missing_value,$variable,d,f, $variable.$yrini-$yrend.$RUN.nc";
system "ncatted -O -a _FillValue,$variable,c,f,-1e30 $variable.$yrini-$yrend.$RUN.nc";   
        

#############
# add $area to file
print "Add $area to file, for weighting\n";
if (! -e "$area$underscore$RUN.nc") {
system "ncks -O -v $area ANN$yrini.$nctag$RUN.nc $area$underscore$RUN.nc";
}
system "ncks -A -a -v $area $area$underscore$RUN.nc $variable.$yrini-$yrend.$RUN.nc";

if ($nctag eq "taij") {
  # taij (sum(mol, C/m2/yr   *  area))*12/1e15
  system "ncap -s '$variable = $variable*$area' $variable.$yrini-$yrend.$RUN.nc dummy2.nc";
  system "ncwa -O -a $lat,$lon -y ttl -d $lat,-90.0,90.0 -d $lon,-180.0,180.0 dummy2.nc dummy3.nc";
  #system "ncap2 -s '$variable = var.total($lat,$lon)' dummy2.nc dummy3.nc";
  #system "ncra -O -y ttl -v $variable dummy2.nc dummy3.nc";
  #system "ncwa -O -a $lat,$lon -d $lat,-90.0,90.0 -d $lon,-180.0,180.0 dummy2.nc dummy3.nc";
  system "ncap -s '$variable = $variable*12/(1e15)' dummy3.nc dummy.nc";
}else{
  system "ncwa -O -a $lat,$lon -w $area -d $lat,-90.0,90.0 -d $lon,-180.0,180.0 $variable.$yrini-$yrend.$RUN.nc dummy.nc";
}

# Remove unwanted terms
system "ncks -O -x -v $lon dummy.nc dummy.nc";
system "ncks -O -x -v $lat dummy.nc dummy.nc";

if (defined($depth)){
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -O -F -d $depth,$ilev,$ilev,1 -v $variable dummy.nc $OutputFileName";
}else{
  system "ncks -O -v $variable dummy.nc $OutputFileName";
}

system "ncrename -O -h -v $variable,glbAvg_ts $OutputFileName";
system "rm -f dummy.nc";


