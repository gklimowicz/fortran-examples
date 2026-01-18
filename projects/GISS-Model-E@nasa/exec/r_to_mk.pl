#!/usr/bin/perl

#OBJ_LIST
#COMPONENTS
#OPTS_giss_LSM

#INPUT_FILES

#RUN_PARAMETERS

#CPP_OPTIONS

#DTFIX ?

# skip garbage at the start of the file (write as comments
while(<>) {
    last if(/Preprocessor *Options/i || /Run *Options/i || /Object *modules/i);
    print "# $_";
}
print "\n";

if ( /Preprocessor *Options/i ) {
    print "CPP_OPTIONS = \n";
    while(<>) {
	last if ( /End *preprocessor *options/i ) ;
	chop;
	s/\#/\\\#/g;
	s/\!/\#/g;
        s/\'/\"/g; # single quotation marks are not supported in CPP section
	print "CPP_OPTIONS += $_\n";
    }
    print "\n";
    # skip till next section
    if (eof) { 
	print STDERR "Error: no 'End preprocessor options'\n"; exit 1;
    }
    while(<>) { last if(/Run *Options/i || /Object *modules/i); }  
}

if ( /Run *Options/i ) {
    print "RUN_OPTIONS = \n";
    while(<>) {
	last if ( /Object *modules/i || /Preprocessor *Options/i ) ;
	chop;
	s/\!/\#/g;
	print "RUN_OPTIONS += $_\n";
    }
    print "\n";
}

#### hack ! hack !
#### repeat parsing of CPP options to allow "Run *Options"
#### to be specified before CPP opts
if ( /Preprocessor *Options/i ) {
    print "CPP_OPTIONS = \n";
    while(<>) {
        last if ( /End *preprocessor *options/i ) ;
        chop;
        s/\#/\\\#/g;
        s/\!/\#/g;
        s/\'/\"/g; # single quotation marks are not supported in CPP section
        print "CPP_OPTIONS += $_\n";
    }
    print "\n";
    # skip till next section
    if (eof) { 
	print STDERR "Error: no 'End preprocessor options'\n"; exit 1;
    }
    while(<>) { last if(/Run *Options/i || /Object *modules/i); }
}

if (eof) { 
    print STDERR "Error: 'Object modules' not found in the rundeck.\n";
    exit 1;
}
if ( /Object *modules/i ) {
    print "OBJ_LIST = \n";
    $OBJ_LIST_O = "";
    while(<>) {
	last if ( /Data input files/i || /Components:/i ) ;
	chop;
	s/\!/\#/g;
	if ( /^\s*$/ ) { print "\n"; next; }
	if ( /^\s*\#/ ) { print "$_\n"; next; }
	$str = $_;
	$str =~ s/\#.*$//;
	$OBJ_LIST_O =join " ", $OBJ_LIST_O, ($str =~ /\w* *\|[^|]*\|/g);
	s/ *\|[^|]*\|//g;  # remove per-file options
	print "OBJ_LIST += $_\n";
    }
    print "\n";
    print "OBJ_LIST_O = $OBJ_LIST_O\n\n"
}

if ( /Components:/i ) {
    print "COMPONENTS = \n";
    while(<>) {
	last if ( /Data input files/i || /Component Options:/i ) ;
	chop;
	s/\!/\#/g;
	s/\s+$//;
	if ( /^\s*$/ ) { print "\n"; next; }
	if ( /^\s*\#/ ) { print "$_\n"; next; }
	print "COMPONENTS += $_\n";
    }
    print "COMPONENTS += profiler\n";
    print "\n";
}

if ( /Component Options:/i ) {
    print "# specific options for components:\n";
    while(<>) {
	last if ( /Data input files/i ) ;
	chop;
	s/\!/\#/g;
	if ( /^\s*$/ ) { print "\n"; next; }
	if ( /^\s*\#/ ) { print "$_\n"; next; }
	print "$_\n";
    }
    print "\n";
}

if (eof) { 
    print STDERR "Error: 'Data input files' not found in the rundeck.\n";
    exit 1;
}
if ( /Data input files/i ) {
    print "INPUT_FILES = \n";
    while(<>) {
	last if ( /Label and Namelist/i) ;
	chop;
	s/\!/\#/g;
	if ( /^\s*$/ ) { print "\n"; next; }
	if ( /^\s*\#/ ) { print "$_\n"; next; }
	print "INPUT_FILES += $_\n";
    }
    print "\n";
    # skip till next section
    while(<>) { last if(/\&\&PARAMETERS/i); }  

}

if (eof) { 
    print STDERR "Error: PARAMETERS not found in the rundeck.\n";
    exit 1;
}
if ( /\&\&PARAMETERS/i ) {
    print "RUN_PARAMETERS = \n";
    while(<>) {
	last if ( /\&\&END_PARAMETERS/i) ;
	chop;
	s/\!/\#/g;
	if ( /^\s*$/ ) { print "\n"; next; }
	if ( /^\s*\#/ ) { print "$_\n"; next; }
	print "RUN_PARAMETERS += $_\n";
    }
    print "\n";
    # skip till next section
    while(<>) { last if(/\&INPUTZ/i); }  

}

# allow skipping INPUTZ if not needed
if (eof) { 
    print STDERR "Warning: INPUTZ not found in the rundeck. Hope it's OK\n";
    exit 0;
}
if ( /\&INPUTZ/i ) {
    print "INPUTZ = \n";
    while(<>) {
	last if ( /\&END/i || /\// ) ;
	chop;
	s/\!/\#/g;
	print "INPUTZ += $_\n";
    }
    print "\n";
}

exit 0 if eof;
# add the rest as comments
while(<>) {
    print "# $_";
}
