#!/usr/bin/perl

while ($_ = $ARGV[0], /^-/) {
  shift;
  last if /^--$/;
  if (/^-par\b/) { $rfile = shift; next;}
  if (/^-in\b/) { $filesource = shift; next;}
  if (/^-out\b/) { $filetarget = shift; next;}
}

print "parfile=$rfile\n";

open(PARFILE, "$rfile") or die "can't open $rfile";

$exist=1;
$lexist=0;
$mtexist=0;
$netcdf=0;

while (<PARFILE>) {
    if ( /End/ ) { last; }
    s/!.*//;
    push @parameter, /([\w.,+_-]+\s*=\s*[\w.,\/+_-]+)/g;
}

foreach $_ ( @parameter ) {
  ($name, $dest) = split /\s*=\s*/;
  if ( $name !~ /^(filedir|remapdir|regridfile|imsource|jmsource|ntilessource|imtarget|jmtarget|ntilestarget|format|intdateline|nfields|title|levels|maintitle)$/ ) {
    print "parameter $name doesn't exist\n";
    $exist=0
    #exit 1;
  } else {
    if ($name =~ "filedir") {
      $filedir=$dest;
    }
    elsif ($name =~ "remapdir") {
      $remapdir=$dest;
    }
    elsif  ($name =~ "regridfile") {
      $regridfile=$dest;
    }
    elsif  ($name =~ "imsource") {
      $imsource=$dest;
    }
    elsif  ($name =~ "jmsource") {
      $jmsource=$dest;
    }
    elsif  ($name =~ "ntilessource") {
      $ntilessource=$dest;
    }
    elsif  ($name =~ "imtarget") {
      $imtarget=$dest;
    }
    elsif  ($name =~ "jmtarget") {
      $jmtarget=$dest;
    }
    elsif  ($name =~ "ntilestarget") {
      $ntilestarget=$dest;
    }
    elsif ($name =~ "format"){
      $format=$dest
    }
    elsif ($name =~ "nfields"){
      $nfields=$dest;
    }
    elsif ($name =~ "title"){
      $title=$dest;
    }
    elsif ($name =~ "levels"){
      $levels=$dest;
      $lexist=1;
    }
    elsif ($name =~ "maintitle"){
      $maintitle=$dest;
      $mtexist=1;
    }
    elsif ($name =~ "intdateline"){
      $idl=$dest;
    }

  }

}

if ( $exist) {

if ( $filesource =~ m/.nc/i) {
   print "netcdf input file \n";
   $netcdf=1;
}
  $files1="$filedir/$filesource";
  print "$files1\n";
  `cp $files1 .`;
  $rfile1="$remapdir/$regridfile";
  print "$rfile1\n";
  `cp $rfile1 .`;

  if ($lexist =~ 0){
    $levels=1;
  }

  if ($mtexist =~ 0){
    $maintitle='no';
  }

if ($netcdf) {
  `./ncll2cs $filesource $filetarget $regridfile $ntilessource $imtarget $jmtarget $ntilestarget $format $idl > output`;
} else {
  `./ll2cs $filesource $filetarget $regridfile $imsource $jmsource $ntilessource $imtarget $jmtarget $ntilestarget $format $nfields $title $levels $maintitle > output`;
}
  open( FILE, "< output" ) or die "Can't open output file : $!";
  while( <FILE> ) {
    print;
  }
  close FILE;
}
