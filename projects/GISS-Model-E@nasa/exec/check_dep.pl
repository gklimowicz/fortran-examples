#!/usr/bin/perl

# this script checks dependencies file for circular dependencies
#
#Usage:  check_dep.pl file_with_dependencies

$filename = shift;

open FILE, $filename or die "check_dep.pl: cannot open $filename";

$str="";
while(<FILE>) {
    chop;
    s/\#.*$//;
    s/^\s*//;
    next if ( ! $_ );
    #print "$_\n";
    $str .= $_;
    if ( $str =~ s/\\.*/ / ) { #"
	 #$str =~ s/\\\s*$"//;
	 #print "oops\n";
	 next;
     }
    #print "str= $str\n";
    if ( $str =~ /^\s*(.+):\s*(.*)/ ) {
	push( @{ $tree{$1}}, split(" ",$2) );
	$str = "";
    }
}

#foreach $k (keys %tree ) {
#    print "key: $k\n";
#    foreach $v ( @{ $tree{$k}} ) {
#	print "--->$v\n";
#    }
#}

print " --- Checking for circular dependencies ---\n";
my $retcode = 0;
foreach $k (keys %tree ) {
    #$k = 'GHY_DRV.o';
    %branch = ();
    if( follow($k) ) {
	print "$k\n";
	print "--->   Circular dependence !!!\n";
	print "\n";
	print "Your source files contain a circular dependence.\n";
	print "The list above shows the dependencies which form a loop.\n";
	print "Please fix this and try to recompile the model.\n";
	print "\n";
	exit 253;
	#retcode = 1;
	#ast;
    }
}

print " ---    dependencies are ok             ---\n";


sub follow{
    my $r = shift;
    my @y;
    my $k;
    #print "root: $r branches: ", join(":", @{ $tree{$r}} ), "\n";

    $branch{$r} = 1;

    @y = @{ $tree{$r}};
    foreach $k (@y) {
	#print "br: $k\n";
	#print "root: $root\n";
	if( $branch{$k} || follow($k) ) {
	    printf "%-33s  - used by:\n", $k; return 1;
	}
    }
    $branch{$r} = 0;
    return 0;
}
