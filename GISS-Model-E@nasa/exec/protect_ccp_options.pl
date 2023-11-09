#!/usr/bin/perl
#
# command to protect CPP options in the rundeck
# when used without options comments out (repalces # with _)
# options between "Preprocessor Options" and "End Preprocessor Options"
#
# when used with -u restores original CPP options
#
# when used with -c also retains the original CPP statement
# (so that it has effect on the rest of the rundeck)
#
# also removes empty lines at the start of rundeck (CPP creates
# extra lines there)

use Getopt::Long;

GetOptions("u","c") || die "problem in GetOptions";

if ($opt_u ) {
    $unprotect = 1;
}

if ($opt_c ) {
    $keep_copy = 1;
}

$start_of_file = 1;
$cpp_protected = 0;

while(<>){

    # skip empy lines at the start of file
    if ($start_of_file && /^\s*$/) { next; }
    $start_of_file = 0;

    # check if we need to protect cpp instructions
    if ( /^\#\#\#cpp_protected/ ) {
	$cpp_protected = 1;
	next;
    }

    if(/^Preprocessor Options/) {
	$inside_cpp_options = 1;
    }
    if(/^End Preprocessor Options/) {
	$inside_cpp_options = 0;
    }

    if( $inside_cpp_options ) {
	if ( (! $cpp_protected) && (! $unprotect) && /^\#/ ) {
	    my $s = $_; 
	    s/^\#/_/;
	    $s =~ s/^(\#\s*define\s+)/$1_/;
	    if ( /^_define\s/ && $keep_copy ) { $_ .= $s; }
	}
	if ( $unprotect ) {
	    if ( /^\s*$/ ) { next; } # remove empty lines lef by #define ...
	    s/^_/\#/;
	}
    }
    print;
}

