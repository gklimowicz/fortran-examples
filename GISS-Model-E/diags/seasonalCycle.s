#!/usr/bin/perl

print "seasonalCycle.s: computes the models global average seasonal cycle at a specified depth\n\n";

###### ---------Script--------- ########
do 'user_input.s';

chdir $DataDir;

$computes = "seasonalCycle_ts";
$new_variablename = "$variable$underscore$computes";
$MonthlyFileName = "$variable.$yrini-$yrend.concattedMonths.$RUN.nc";
$OutputFileName = "$variable.$yrini-$yrend.$computes.lev$ilev.$RUN.nc";
print "Output Name $OutputFileName\n";

# Create month.nc
if (! -e month.nc){
system "ncks -v mon $ObsDir$ObsFilename month.nc";
}

@months = qw(JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC);
print "Averaging annual cycle for years $yrini to $yrend \n";

foreach $month (@months) 
{
   system "scaleacc $month$yrini-$yrend.acc$RUN.nc $nctag";
   system "ncks -O -v $variable $month$yrini-$yrend.$nctag$RUN.nc dummy_$month.nc";
}
system "ncecat dummy_{JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC}.nc $MonthlyFileName";
system "rm -f dummy_*.nc";

#Rename missing_value to FillValue before computation
system "ncatted -O -a missing_value,$variable,d,f, $MonthlyFileName";
system "ncatted -O -a _FillValue,$variable,c,f,-1e30 $MonthlyFileName";   

#############
if ($area eq "oxyp") {
system "scaleacc JAN$yrini-$yrend.acc$RUN.nc oij";
system "ncks -O -v $area JAN$yrini-$yrend.oij$RUN.nc $area$underscore$RUN.nc";
}else{
system "scaleacc JAN$yrini-$yrend.acc$RUN.nc aij";
system "ncks -O -v $area JAN$yrini-$yrend.aij$RUN.nc $area$underscore$RUN.nc";
}

if (defined($depth)){
  # extract $depth at srf
  print "Extracting depth at level $ilev \n";
  system "ncks -O -F -d $depth,$ilev,$ilev,1 -v $variable $MonthlyFileName dummy.nc";
} else  {
  system "ncks -O -v $variable $MonthlyFileName dummy.nc";
}

system "ncks -A -a -v $area $area$underscore$RUN.nc dummy.nc";
system "ncwa -h -O -w $area -a $lat,$lon dummy.nc $OutputFileName";

# Remove unwanted terms
system "ncks -O -x -v $lon $OutputFileName $OutputFileName";
system "ncks -O -x -v $lat $OutputFileName $OutputFileName";
system "ncks -O -x -v $area $OutputFileName $OutputFileName";
system "ncatted -O -a long_name,$variable,d,s, $OutputFileName";

system "ncrename -O -h -v $variable,$new_variablename $OutputFileName";
if (defined($depth)){
system "ncrename -O -h -v zoc,dep -d zoc,dep $OutputFileName";
}
system "ncrename -d record,mon $OutputFileName";
system "ncks -A -v mon month.nc $OutputFileName";
system "ncatted -O -a _FillValue,$new_variablename,d,f, $OutputFileName";

# Computes the northern and southern hermisphere averages
if ($variable eq "oicefr"){
   $NH="NortH";
   $SH="SoutH";
   system "ncwa -h -O -w $area -a $lat,$lon --mask_condition '$lat>0' dummy.nc dummy1.nc";
   system "ncrename -v $variable,$new_variablename$underscore$NH dummy1.nc";
   $OutputFileName = "$variable.$yrini-$yrend.$computes.lev$ilev.$RUN.$NH.nc";
   system "ncks -h dummy1.nc $OutputFileName";
   system "ncwa -h -O -w $area -a $lat,$lon --mask_condition '$lat<0' dummy.nc dummy2.nc";
   $OutputFileName = "$variable.$yrini-$yrend.$computes.lev$ilev.$RUN.$SH.nc";
   system "ncrename -v $variable,$new_variablename$underscore$SH dummy2.nc";
   system "ncks -h dummy2.nc $OutputFileName";
   print "$OutputFileName \n";
}

system "rm -f dummy.nc";



