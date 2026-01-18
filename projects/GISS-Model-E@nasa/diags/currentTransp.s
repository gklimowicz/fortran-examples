#!/usr/bin/perl

# Computes the transport of ACC, Gulf Stream and Kuroshio

###### ---------Script--------- ########

do 'user_input.s';
print "\n";

chdir $DataDir;
#$OutputFileName = "$variable.$yrini-$yrend.glbAvg_lev$ilev.$RUN.nc";
#print "Output Name $OutputFileName\n";

system "echo Current Transports for $RUN from $yrini to $yrend > dummy.txt";
print " \n";
$year = $yrini;
while ($year <= $yrend) 
{
  system "scaleacc ANN$year.acc$RUN.nc oij,ojl";
  system "prtostat $RUN ANN$year >> dummy.txt ";
  system "grep 'Gulf Stream' dummy.txt >  GulfStream$yrini-$yrend.txt";
  system "grep 'Kuroshio' dummy.txt >  Kuroshio$yrini-$yrend.txt";
  system "grep 'Drakes' dummy.txt >  DrakesPassage$yrini-$yrend.txt";
  $year=$year + 1;
}


