#!/usr/bin/perl

# adapt the rundeck to local environment

# get modelErc
if (exists $ENV{MODELERC}) {
    $modelerc = $ENV{MODELERC}
}
else {
    $modelerc = (getpwuid($>))[7]."/.modelErc";
}

#defaults
$PNETCDFHOME = "";

if ( -f $modelerc ) {
    open MODELERC, $modelerc or die "can't open $modelerc";
    while(<MODELERC>) {
	$PNETCDFHOME = $1 if /^ *PNETCDFHOME *= *(\S+)/;
    }
    close MODELERC;
}

while(<>) {
    if ( (! $PNETCDFHOME) && /^\s*OPTS_dd2d\s*=[^#]*NC_IO=PNETCDF/  ) { #/
	 print "!make> no PNETCDFHOME - removing NC_IO=PNETCDF\n";
	 # get rid of PNETCDFHOME
	 s/NC_IO=PNETCDF//g;
	 # if no other options -remove the line
	 s/OPTS_dd2d\s*=\s*(#.*)*//;
    }
    print ;

}
