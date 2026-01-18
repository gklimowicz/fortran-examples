#!/usr/bin/perl -w
#
# HELP: get_out4archive.pl -h
# SCRIPT: get_out4archive.pl           
# get modelE acc or hycom out files from archive
# AUTHOR: Nick Tausnev, ntausnev@giss.nasa.gov 
#
# PURPOSE: copy modelE  *.acc$RunId OR *.out$RunId files 
# from  directories ${ARCHIVE}/RunId/00ACC OR ${ARCHIVE}/RunId/00HYC
#    to 
# user setting  $outDir = "..."; if $outDir is not in arguments
#  DEFAULT:  $outDir = "/discover/nobackup/ntausnev/RUNS_ME/${RunId}/00ACC";
#            $outDir = "/discover/nobackup/ntausnev/RUNS_ME/${RunId}/00HYC";
#
# Problems : not tested if some archive files absent?
# SOME hard coded  variables user should adjust
#          
# Script assumes that archive data 
# ${ARCHIVE}/${RunId}/${ARCH_SUB_DIR}
# ${ARCHIVE}/${RunId}/00ACC or ${ARCHIVE}/${RunId}/00HYC
# Data in form: 
# ${Field}${RunId}_\?\?\?\?-\?\?\?\?\.tar ----> accEh127f9c_1850-1859.tar or
#                                               outEh127f9c_1850-1859.tar
#
#
# Need some user updating 
#
#
use warnings;
use strict;
use File::Basename;
use File::Copy;
use Getopt::Long;
#
use Env qw($ARCHIVE);
#
my $command = $0;
$command =~ s#^.*/##;

my $outDir;

my $help  = '';
my $RunId = '';
my $Field = '';
my $dir_out = '';

my $startYear = -1;
my $lastYear  = -1;

#call the function passing the param name and variable to assign the value to

GetOptions ('h|help'        => \$help,
            'r|run=s'       => \$RunId,
            'f|field=s'     => \$Field,
            's|startYear=i' => \$startYear,
            'l|lastYear=i'  => \$lastYear,
            'd|dir_out=s'   => \$dir_out);
#
##########################################################
########### DEFINE FILES AND VARIABLES HERE ##############
##########################################################
#
if ($help){
print "The help parameter was supplied to the program\n";
print <<"EOT";
Script can be used with next options:
  $command -h
  $command -r RunId -f acc|out -s startYear -l lastYear  [ -d dir_out ]

EOT
exit (0);
}

if ( ! $RunId ){
  print "RunId was not supplied\n";
  exit (1)
} else {
  print "RunId set to $RunId\n";
}

if ( ! $Field ){
  print "Field was not supplied: avilable  acc OR out \n";
  exit (1)
} else {
  print "Field set to $Field\n";
}

if ( $startYear == -1 ){
  print "Start year was not supplied\n";
  print "\$startYear = $startYear \n";
} else {
  print "Start year set to $startYear\n";
  if ( $startYear  < 0  ) {
     print "\nStart year  should be >= 0\n";
     print "!!! Exit from script $command \n";
     exit(1);
  }
}

if ( $lastYear == -1 ){
  print "Last year was not supplied\n";
  $lastYear = $startYear;
} else {
  print "Last  year set to $lastYear\n";
  if ( $lastYear  < 0  ) {
     print "\nLast  year  should be >= 0\n";
     print "!!! Exit from script $command \n";
     exit(1);
  }
  my $diff = $lastYear - $startYear;
  if ( $diff < 0  ) {
     print "\nLast year ( $lastYear ) should be >= Start year ( $startYear )\n";
     print "!!! Exit from script $command \n";
     exit(1);
  }
}

if ( ! $dir_out ){
  $outDir = "/discover/nobackup/ntausnev/RUNS_ME/${RunId}/z_00HYC";    ### DEFAULT
  print "outDir was not supplied it will use default value:\n  $outDir\n";
} else {
  $outDir = $dir_out;
  print "outDir set to $dir_out\n";
}

my $ARCH_SUB_DIR = "";

if ( $Field eq "acc" ) {
     $ARCH_SUB_DIR = "00ACC";
}  
if ( $Field eq "out" ) {
     $ARCH_SUB_DIR = "00HYC";
}  
my $archiveDir = "${ARCHIVE}/${RunId}/${ARCH_SUB_DIR}";  ### Hard coded !!! 
print "Files will be copied from remote dir=${archiveDir}\n";

opendir (Archive, $archiveDir) or die "Can't open directory=$archiveDir: $!";
my @files = grep { /${Field}${RunId}_\?\?\?\?-\?\?\?\?\.tar$/ } readdir(Archive);
my @files1 = grep { /\.tar$/ } readdir(Archive);
close(Archive);

my $file;
my $locPathFile;
my ($year1, $year2);
my @files2copy;
my $flag =0;    # if $flag = 1  add name of file at @files2copy   
@files = glob("$archiveDir/${Field}${RunId}_\?\?\?\?-\?\?\?\?\.tar");  ### Hard coded (out)
foreach $file (@files) {
# print "$file\n";
  $locPathFile = basename($file);
# print "locPathFile=$locPathFile \n";
  $locPathFile = $locPathFile;
  if( $locPathFile =~ /(\d\d\d\d)-(\d\d\d\d)/){
     $year1 = "$1";
     $year2 = "$2";
#    print "year1,year2 = $year1,$year2 \n";  # comment !!!
    }
    else{
     print "Invalid years at the name ${Field} files!\n";
  }  
  if ( ($startYear - $year1 >= 0) && ( $year2 - $startYear >= 0)  ){
     print " startYear=$startYear is between $year1 and $year2 \n";
     $flag = 1;
  }
  if ( $flag == 1 ) {
     push( @files2copy, $file);
  }  
  if ( ($lastYear - $year1 >= 0) && ( $year2 - $lastYear >= 0)  ){
     print " lastYear=$lastYear is between $year1 and $year2 \n";
     last;
  }
}
my $dir;
my $prefix;

chdir ${outDir} or die "Can't cd to ${outDir} $!\n";

print "\nThe list of files to copy from archive:\n";
print map {"$_ \n" } @files2copy;
foreach $file (@files2copy) {
  print "file=$file\n";
  $locPathFile = basename($file);
  ($dir,$prefix) = split(/\./,$locPathFile);
  
  copy($file, ${locPathFile}) or die "copy failed: $!";
  my $tar = `tar xvf ${locPathFile} 2>&1`;
  my $remove_tar = `  rm -f ${locPathFile}`;
  print "\nUntar data ${locPathFile}\n";
# my $moveUp = ` mv $dir/* . ; rm -f ${locPathFile}`;
}

# End of script
exit (0);

