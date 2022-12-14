#!/usr/bin/perl

#### so that there is not confusion over the years averaged

use Cwd;
$myDir = getcwd;
print "current directory $myDir \n";

##### -------- INPUTS ------- #####
chdir $myDir;
do 'user_input.s';
print "\n";

##### -------- Script ------- #####
##### -------- averaging files, annual means, climatologies, seasonal means  ------- #####
# copies over accfiles, sumfiles for annual means, climaological means
# and seasonal means, seasonal cycles and then deletes accfiles

=pod
chdir $myDir;
print "Doing avgACC.s \n";
do 'avgACC.s';
print "\n";
=cut

##### -------- mean annual cycle at a certain level  ------- #####
chdir $myDir;
print "nctag is $nctag\n";
print "Doing seasonalCycle.s \n";
do 'seasonalCycle.s';
print "\n";
$FileName1 = "$variable.$yrini-$yrend.$computes.lev$ilev.$RUN";

##### -------- mean annual cycle from obs    ------- #####
chdir $myDir;
print "Doing seasonalCycleObs.s \n";
print "$ObsDir$ObsFilename \n";
do 'seasonalCycleObs.s';
print "\n";
$new_variablename1 = "$new_variablename";
$FileName2 = "$variable_obs.$computes.lev$ilev.$ObsFilename";
substr($FileName2,rindex $FileName2,'.') = '';

##### -------- mean annual cycle diff from obs    ------- #####
chdir $myDir;
print "Doing diff_seasons.s \n";
do 'diff_seasons.s';
print "\n";

print "datadir =  $DataDir \n";
print "new_variablename = $new_variablename1 \n";
print "filename1 =  $FileName1 \n";
print "filename2 =  $FileName2 \n";

###### ---------PYTHON Script--------- ########
##invoke the python script to plot model seasonal cycle
chdir $myDir;

if ($variable ne "oicefr"){
    system "python3 plot_line.py $new_variablename1 $DataDir$FileName1.nc $DataDir$FileName2.nc";

} else{
    system "python3 plot_line.py $new_variablename1$underscore$NH $DataDir$FileName1.$NH.nc $DataDir$FileName2.$NH.nc";

    system "python3 plot_line.py $new_variablename1$underscore$SH $DataDir$FileName1.$SH.nc $DataDir$FileName2.$SH.nc";
} 

###### ------------------------------- ########
##### -------- climatology maps at certain level  ------- #####
chdir $myDir;
print "Doing clim_map.s \n";
do 'clim_map.s';
print "\n";
$OutputFileName1 = "$OutputFileName1";

##### -------- maps from observations   ------- #####
chdir $myDir;
print "Doing obs_map.s \n";
do 'obs_map.s';
print "\n";
$OutputFileName2 = "$OutputFileName";

chdir $myDir;
print "Doing diff_maps.s \n";
do 'diff_maps.s';
print "\n";
$OutputFileName3 = "$OutputFileName";

print "variable = $variable \n";
print "outputfile1 = $DataDir$OutputFileName1 \n";
print "variable_obs = $variable_obs \n";
print "outputfile2 = $DataDir$OutputFileName2 \n";
print "variable_diff = $variable$underscore$diff  \n";
print "outputfile3 = $DataDir$OutputFileName3 \n";
print "model run = $RUN \n";

###### ---------PYTHON Script--------- ########
##invoke the python script
chdir $myDir;

system "python3 plot_map.py $variable $DataDir$OutputFileName1 $DataDir$OutputFileName2 $DataDir$OutputFileName3";
###### ------------------------------- ########

##### -------- global averaged timeseries at a certain level   ------- #####
chdir $myDir;
print "Doing map_glbavg_ts.s \n";
do 'map_glbavg_ts.s';
print "\n";

###### ---------PYTHON Script--------- ########
##invoke the python script to plot model seasonal cycle
chdir $myDir;
system "python3 plot_linets.py glbAvg_ts $DataDir$OutputFileName";
###### ------------------------------- ########

##### -------- vertical sections in different basins   ------- #####
##### -------- diff basin vertical sections from obs   ------- #####
chdir $myDir;
print "Doing basinAvg_obs.s \n";
do 'basinAvg_obs.s';
print "\n";
$OutputFileName_obs="$DataDir$OutputFileName";

chdir $myDir;
print "Doing basinAvg_model.s \n";
do 'basinAvg_model.s';
print "\n";
$OutputFileName_model="$DataDir$OutputFileName";

chdir $myDir;
print "Doing diffBasinAvg.s \n";
do 'diffBasinAvg.s';
print "\n";

###### ---------PYTHON Script--------- ########
##invoke the python script
chdir $myDir;
if ($area eq "oxyp"){
system "python3 plot_section.py $variable $DataDir$OutputFileName";
}else{
$ann="ann";
print "variable $variable \n";
print "variable obs $variable_obs \n";
print "file model $OutputFileName_model \n";
print "file obs $OutputFileName_obs \n";

system "python3 plot_surf_basin_avg.py $variable $OutputFileName_model $variable_obs $OutputFileName_obs";
}
###### ------------------------------- ########

##### -------- MOC vertical sections, all basins      ------- #####
  chdir $myDir;
  print "Doing basinAvg_MOC.s \n";
  do 'basinAvg_MOC.s';
  print "\n";

###### ---------PYTHON Script--------- ########
##invoke the python script
chdir $myDir;
print "$DataDir$OutputFileName";
system "python3 plot_moc.py $DataDir$OutputFileName";
###### ------------------------------- ########

##### -------- AMOC max at 26N timeseries              ------- #####
  chdir $myDir;
  print "Doing AMOCvertMax.s \n";
  do 'AMOCvertMax.s';
  print "\n";

###### ---------PYTHON Script--------- ########
##invoke the python script to plot model seasonal cycle
chdir $myDir;
system "python3 plot_linets.py sf_Atl $DataDir$OutputFileName";

###### ------------------------------- ########


##### -------- Current Transport timeseries ------- #####
## transports for Kuroshio, Gulf Stream and ACC
  chdir $myDir;
  print "Doing currentTransp.s \n";
  do 'currentTransp.s';
  print "\n";

###### ---------PYTHON Script--------- ########
##invoke the python script to plot model seasonal cycle

chdir $myDir;
$GS = GulfStream;
$KS = Kuroshio;
$DP = DrakesPassage;
system "python3 plot_text.py $DataDir$GS$yrini-$yrend.txt $DataDir$KS$yrini-$yrend.txt $DataDir$DP$yrini-$yrend.txt";
###### ------------------------------- ########
