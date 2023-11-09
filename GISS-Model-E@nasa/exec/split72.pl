#!/usr/bin/perl
#@auth I. Aleinov
# This script splits lines of fortran code so that no line is
# longer than 72 characters
# Bugs: splits character strings in the middle, which may cause
# problems with some compilers and in general is not a very 
# good idea.

#for debugging:
$num_long_lines = 0;
$num_split_lines = 0;

while(<>) {

    # if comments or line < 73 characters do nothing
    if ( /^[Cc]/ || /^\s*!/ || ! /.{73}/ ) { print; next; }

    chop;

    # find comment at the end of line (if any)
    $comment = $_;
    $counter1 = 0;
    while( $comment ) {
	$comment =~ s/^[^!^\"^\']+//;
	$comment =~ s/^\".*?\"//;
	$comment =~ s/^\'.*?\'//;
	if ( $comment =~ /^!/ ) { last; }
	$counter1 ++;
	if ( $counter1 > 64 ) { #assume no more than 32 char strings / line
	    print STDERR "split72: unmatched quotation mark\n";
	    $comment = ""; # assume no comment
	    last;
	    }
	}

    # remove comment from the end of the line
    if ( $comment ) {
	$comment_length = length($comment);
	s/.{$comment_length}$//;
    }

    #print "string:\n{$_}X\ncomment:\n{$comment}X\n";
    
    $counter = 0;
    while ( /(.{73})/ ) {
        s/^(.{1,72})\b(?![%])/     &/;  #"# don't split on chars in []
        print "$1\n";
        $counter++;
        if ( $counter > 32 ) {
            print STDERR "split72: Can't split the line or line too long\n";
            exit (1);
        }
    }
    print "$_$comment\n";

    $num_long_lines ++;
    if ( $counter ) { $num_split_lines ++; }

}

print STDERR "split72: Info:\n";
print STDERR "    lines longer than 72 chars: $num_long_lines\n";
print STDERR "    lines split $num_split_lines\n";
