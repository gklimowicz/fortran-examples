#!/usr/bin/perl

# creates a 2D plot difference of model - gridded observations
# from a specified level

###### ---------Script--------- ########
print "Doing diff_maps. Subtract model 2D climatology map from 2D obs map.\n";

do 'user_input.s';

if (index($variable_obs, "_ann") != -1) {
     print "'$variable_obs' contains ann.\n";
     $variable_obs =~ s/ann/mon/g;   
}

chdir $DataDir;
$diff="diff";
$new_variable = "$variable$underscore$diff";
print "$new_variable \n";

@data_array = ("DJF","JJA","ANN");
for my $months (@data_array){
  $yrin = $yrini;
  if ($months eq "DJF"){
    $yrin = $yrini+1;
  }
  if ($months eq "ANN"){
     if (index($variable_obs, "_mon") != -1) {
     print "'$variable_obs' contains mon.\n";
     $variable_obs =~ s/mon/ann/g;   
  }
  }
  $FirstFileName = "$variable.$months$yrin-$yrend.map_lev$ilev.$RUN.nc";
  print "First file $FirstFileName \n";
  $SecondFileName = "$variable_obs.$months.map_lev$ilev.$RUN2.nc";
  print "Second file $SecondFileName \n";

  $OutputFileName = "$new_variable.$months$yrin-$yrend.diffMap_lev$ilev.$RUN$underscore$RUN2.nc";
   print "output file $OutputFileName \n";
 
  # Make all variables have the same name to use ncdiff
  #system "ncrename -v $variable_obs,$new_variable $SecondFileName dummy2.nc";  
  system "ncrename -v $variable,$new_variable $SecondFileName dummy2.nc";  
  system "ncrename -v $variable,$new_variable $FirstFileName dummy1.nc";                      
=pod
  #rename lat/lon from lato/lono
  if ($lat_obs ne "lat") {
    system "ncrename -O -v $lat_obs,lat -d $lat_obs,lat dummy2.nc";  
    system "ncrename -O -v $lon_obs,lon -d $lon_obs,lon dummy2.nc";
  }
  if ($lat ne "lat") {
    system "ncrename -O -v $lat,lat -d $lat,lat dummy1.nc";  
    system "ncrename -O -v $lon,lon -d $lon,lon dummy1.nc";  
  }
=cut

  system "ncdiff -v $new_variable dummy1.nc dummy2.nc $OutputFileName";
  system "rm dummy*";
}

