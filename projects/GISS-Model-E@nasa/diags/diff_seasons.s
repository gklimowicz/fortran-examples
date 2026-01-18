#!/usr/bin/perl

# creates a 2D plot difference of model - gridded observations
# from a specified level

do 'user_input.s';

$diff = "diff";
$computes = "seasonalCycle_ts";
$old_variablename = "$variable$underscore$computes";
$new_variablename = "$old_variablename$underscore$diff";

$FirstFileName = "$variable.$yrini-$yrend.$computes.lev$ilev.$RUN.nc";
print "First file $FirstFileName \n";

$SecondFileName = "$variable_obs.$computes.lev$ilev.$RUN2.nc";
print "Second file $SecondFileName \n";

###### ---------Script--------- ########
chdir $DataDir;

$OutputFileName = "$variable.$yrini-$yrend.$computes.Diff_lev$ilev$RUN$underscore$RUN2.nc";

# Make all variables have the same name to use ncdiff
system "ncrename -v $old_variablename,$new_variablename $SecondFileName dummy2.nc";  
system "ncrename -v $old_variablename,$new_variablename $FirstFileName dummy1.nc";                      

# Add mon variable to model
system "ncks -A -a -v mon dummy2.nc dummy1.nc";

system "ncdiff -v $new_variablename dummy1.nc dummy2.nc $OutputFileName";
system "rm dummy*";

