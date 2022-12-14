#!/usr/bin/perl 
#
# parse .depend file for each component and create a file with
# dependencies between the components
#

$verbose = 0;

%mod_provided = ();
%mod_needed = ();
%comp_depends = ();

if ( $#ARGV < 0 ) {
    #print_main_index();
    exit;
}


while( $comp = shift ) {
    my $file_to_open = $comp."/.depend";
    $comp .= "_dir";
    open(SRCFILE, $file_to_open) or die "can't open $file_to_open\n";
    print "parsing $file_to_open\n";
    parse_file();
    close SRCFILE; #just in case want to reset line counter ($.)
}

foreach $comp ( keys %mod_needed ) {
    foreach $mod_name ( keys %{$mod_needed{$comp}} ) {
	$comp_depends{$comp}{ $mod_provided{$mod_name} } = 1;
	print "MOD NEEDED $comp -- $mod_name, : $mod_provided{$mod_name}\n" 
	    if $verbose;
    }
}

open ( DEP, ">.depend_subdirs" ) || die "can't open $filename";
foreach $comp ( keys %comp_depends ) {
    print DEP "$comp: ";
    foreach $depends_on ( keys %{$comp_depends{$comp}} ) {
	print DEP " $depends_on" if $depends_on ne $comp;
    }
    print DEP "\n";
}
	


sub parse_file {
    while( <SRCFILE> ) {
	chop;
	if ( /\s*(\w+\.mod): +(\S+\.smod)/i ) { $mod_provided{$1} = $comp; next; }
	if ( /\s*(\w+\.o): +/ ) {
	    my $obj_file = $1;
	    while ( s/(\w+\.mod)//i ) { $mod_needed{$comp}{$1} = 1; }
	}

    }
}

