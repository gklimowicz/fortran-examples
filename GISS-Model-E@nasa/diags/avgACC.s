#!/usr/bin/perl

# copies over accfiles, sumfiles for annual means, climaological means
# and seasonal means, then deletes accfiles

###### ------------------------------------------ ########
###### --------- copy acc files --------- ########
do 'user_input.s';

print "Copying acc files, $RUN years from $yrini to $yrend \n";
print " \n";

chdir $DataDir;

$year = $yrini;
while ($year <= $yrend)
{
  print "Copying year $year.\n";
  system "cp $ACCDataDir/*$year.acc$RUN.nc .";
  $year=$year + 1;
}


###### ------------------------------------------ ########
###### --------- annual mean files --------- ########
print "ANN mean: averaging over each year, $RUN years from $yrini to $yrend \n";
print " \n";

chdir $DataDir;

$year = $yrini;
while ($year <= $yrend) 
{
  print "Doing year $year.\n";
  system "sumfiles {JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC}$year.acc$RUN.nc";
  $year=$year + 1;
}

###### ------------------------------------------ ########
###### --------- climatological means --------- ########
print "ANNyrini-yrend: averaging over each long term climatology, $RUN from $yrini to $yrend \n";
print " \n";

@filenames = ("ANN$yrini.acc$RUN.nc");
 #print $_, "\n" for @filenames;
 $year = $yrini+1;
 while ($year <= $yrend)
 {
 push @filenames, "ANN$year.acc$RUN.nc";
 $year=$year + 1;
 }
 print $_, "\n" for @filenames;
 system "sumfiles @filenames";

###### ------------------------------------------ ########
###### --------- seasonal climatologies --------- ########
print "DJF and JJA: averaging seasonal climatologies, $RUN from $yrini to $yrend \n";
print " \n";

system "sumfiles {DEC,JAN,FEB}*.acc$RUN.nc";
system "sumfiles {JUN,JUL,AUG}*.acc$RUN.nc";

###### ------------------------------------------ ########
###### --------- annual cycles          --------- ########
@months = qw(JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC);
print "Annual Cycle: averaging mean annual cycle for years, $RUN $yrini to $yrend \n";
print " \n";

foreach $month (@months)
{
 print "Doing $month \n";
 @filenames = ("$month$yrini.acc$RUN.nc");
 #print $_, "\n" for @filenames;
 $year = $yrini+1;
 while ($year <= $yrend)
 {
 push @filenames, "$month$year.acc$RUN.nc";
 $year=$year + 1;
 }
 print $_, "\n" for @filenames;
 system "sumfiles @filenames";
}

###### ------------------------------------------ ########
###### --------- remove original accfiles-------- ########
# remove .acc files
#system "rm -R -f {JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC}*.acc$RUN.nc";
