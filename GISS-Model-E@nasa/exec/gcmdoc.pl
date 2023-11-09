#!/usr/bin/perl
#This script parses fortran files and creates documentation
#in HTML format. It utilizes information provided with @.* tags.
#Author: I. Aleinov
#Version:  1.0
#Usage:  gcmdoc.pl [-D output_dir] file1.f [file2.f ...]
#           output_dir - a directory were to write HTML files

use Getopt::Long;         #module to process command line options
use Shell qw(date cp ls);    #module to use shell commands

# defaults for some global definitions:

$output_dir = ".";
$doc_dir = "../doc";
$run_name = "";
$component_name = "";

#some useful patterns

# this pattern is not very precise, but if one makes it absolutely
# correct it will be a mile long with this crazy fortran syntax
$some_decl = "integer|real|logical|complex";
$some_decl = "(?:$some_decl)(?:\\s*\\*\\s*(?:\\d+|\\(\\[0-9 ]+\\)))?";
$some_decl .= "|double\\s+precision";
#$some_decl .= "|character(?:\\s*\\*\\s*(?:\\d+|\\([0-9a-z*= ]+\\)))?";
$some_decl .= "|character(?:\\s*\\*?\\s*(?:\\d+|\\([0-9a-z*= ]+\\)))?";
$some_decl .= "|type\\s*\\(\\s*\\w+\\s*\\)";
$some_decl .= "|use\\s*\\w+\\s*,\\s*only\\s*:";
# $some_decl .= "|interface\\s+\\w+";

$parenth  = "\\([^()]*\\)";
$parenth2 = "\\((?:[^()]|$parenth)*\\)";
$parenth3 = "\\((?:[^()]|$parenth2)*\\)";

#print "SOME DECL >>>$some_decl<<<\n";

#an example
#GetOptions("s", "e=s", "f=s", "I=s@", "m=s", "c", "p", "g", "h", "o=s", "a=s")
#                       || die "problem in GetOptions";

GetOptions("D=s", "O=s", "R=s", "CPP=s", "C=s") || die "problem in GetOptions";

if( $opt_O ) {
    $output_dir = $opt_O;
    $output_dir =~ s/\/?\s*$//;
}

if( $opt_D ) {
    $doc_dir = $opt_D;
    $doc_dir =~ s/\/?\s*$//;
}

if( $opt_R ) {
    $run_name = $opt_R;
}

if( $opt_C ) {
    $component_name = $opt_C;
}


print "will write output to: $output_dir\n";

# define some useful global variablers

$current_date = date();
$pwd = `pwd`; chop $pwd;
$abs_doc_dir = $doc_dir; $abs_doc_dir =~ s/^(?!\/)/$pwd\//;
$abs_output_dir = $output_dir; $abs_output_dir =~ s/^(?!\/)/$pwd\//;
$path_doc_output = relpath($abs_doc_dir,$abs_output_dir);
$path_output_doc = relpath($abs_output_dir,$abs_doc_dir);
if ( ! $run_name ) { $run_name = $output_dir; }

# the following is used to generate "official" web pages
$html_include_dir = $ENV{'HTMLINCLUDEDIR'};

#init global hashes:
%db_vars = ();       #var name: module:sub:var
%db_subs = ();       #sub name: module:sub
%db_modules = ();
%db_files = ();


if ( $#ARGV < 0 ) {
    print_main_index();
    exit;
}


while( $current_file = shift ) {
    my $file_to_open = $current_file;
    #if -cpp option is specified try to open .cpp instead of .f
    if ( $opt_CPP ) { $file_to_open =~ s/$/$opt_CPP/; }
    open(SRCFILE, $file_to_open) or die "can't open $file_to_open\n";
    print "parsing $file_to_open\n";
    parse_file();
    close SRCFILE; #just in case want to reset line counter ($.)
}

postprocess_lists();

# printing output to html files

print  "printing file list\n";
htm_start("$output_dir/files.html","GCM Source Files");
print HTM '<H3><Center>FILES</Center></font></H3>'."\n";
print HTM "<dl>\n";
foreach $name ( sort keys %db_files ) {
    print HTM "<dt><b>";
    htm_link($name,"$name.html");
    print HTM "</b><BR>\n";
    #print "printing $name\n";
    #print HTM "$name\n";
    #print HTM "<UL>\n";
    #print HTM "Summary: $db_files{$name}{sum}<BR>\n";
    print HTM "<dd>Modules: ";
    #while ( $name_mod = pop @{$db_files{$name}{modules}} ) {
    foreach $name_mod ( @{$db_files{$name}{modules}} ) {
	#print HTM " $name_mod";
	print HTM " ";
	htm_link("$name_mod","$name_mod.html");
    }
    print HTM "<BR>\n";
    #print HTM "Summary: $db_files{$name}{sum}<BR>\n";
    htm_text("Summary: $db_files{$name}{sum}");
}
print HTM "</dl>\n";
htm_end();

print  "printing module list\n";
htm_start("$output_dir/modules.html","GCM Fortran Modules");
print HTM '<H3><Center>MODULES</Center></font></H3>'."\n";
print HTM "<dl>\n";
foreach $name ( sort keys %db_modules ) {
    print HTM "<dt><b>";
    htm_link($name,"$name.html");
    print HTM "</b><BR>\n";
    #print HTM "<dd>$db_modules{$name}{sum}<BR>\n";
    print HTM "<dd>";
    htm_text("$db_modules{$name}{sum}");
}
print HTM "</dl>\n";
htm_end();

print  "printing variables\n";

htm_start("$output_dir/vars_db.html","GCM ALL Variables");
print HTM '<H3><Center>List of Database Variables</Center></font></H3>'."\n";
@var_list = ();
foreach $name ( keys %db_vars ) {
    if ( $db_vars{$name}{tag} =~ /\bdbparam\b/i ) {
	push @var_list, $name;
    }
}
htm_prt_vars( "sl", @var_list );
htm_end();

htm_start("$output_dir/vars_nl.html","GCM ALL Variables");
print HTM '<H3><Center>List of Namelist Variables</Center></font></H3>'."\n";
@var_list = ();
foreach $name ( keys %db_vars ) {
    if ( $db_vars{$name}{tag} =~ /\bnlparam\b/i ) {
	push @var_list, $name;
    }
}
htm_prt_vars( "sl", @var_list );
htm_end();

htm_start("$output_dir/vars_global.html","GCM ALL Variables");
print HTM '<H3><Center>List of All Global Variables</Center></font></H3>'."\n";
@var_list = ();
foreach $name ( keys %db_vars ) {
    #if ( $db_vars{$name}{tag} =~ /\bnlparam\b/i ) {
    if ( $name =~ /\w+::\w+/ ) {
	push @var_list, $name;
    }
}
htm_prt_vars( "sl", @var_list );
htm_end();


htm_start("$output_dir/vars_all.html","GCM ALL Variables");
print HTM '<H3><Center>List of All Variables</Center></font></H3>'."\n";
htm_prt_vars( "sl", keys %db_vars );
htm_end();

print  "printing files\n";
# create USE dependencies list
@list_db_files = sort keys %db_files;
foreach $name ( keys %db_files ) {
    foreach $mod_name ( keys %{$db_files{$name}{depend_on_mod}} ) {
	$db_files{$name}{depend_on_file}{$db_modules{$mod_name}{file}} = 1;
    }
}


foreach $name ( keys %db_files ) {
    htm_start("$output_dir/$name.html","$name");
    print HTM "<H3><Center>$name</Center></font></H3>\n";
    #$str = $db_files{$name}{sum};
    #$str =~ s/>/&gt;/g; $str =~ s/</&lt;/g; 
    #print HTM "Summary: $str<BR>\n";
    htm_text("Summary: $db_files{$name}{sum}");
    print HTM "Author : $db_files{$name}{auth}<BR>\n";
    print HTM "Version: $db_files{$name}{ver}<BR>\n";
    if ( $db_files{$name}{usage} ) { htm_usage($db_files{$name}{usage}); }
    print HTM "<HR width=10%>\n";
    print HTM "Modules: \n";
    print HTM "<dl>\n";
    foreach $mod_name ( @{$db_files{$name}{modules}} ) {
	print HTM "<dt>";
	htm_link($mod_name,"$mod_name.html");
	print HTM "<br>";
	#print HTM "<dd>$db_modules{$mod_name}{sum}<BR>\n";
	print HTM "<dd>";
	htm_text("$db_modules{$mod_name}{sum}");
    }
    print HTM "</dl>\n";
    print HTM "<HR width=10%>\n";
    print HTM "Global Subroutines: \n";
    print HTM "<dl>\n";
    foreach $sub_name ( sort @{$db_files{$name}{subs}} ) {
	print HTM "<dt>";
	htm_link($sub_name,":$sub_name.html");
	print HTM "<br>";
	$sub_name = ":$sub_name";
	#print HTM "<dd>$db_subs{$sub_name}{sum}<BR>\n";
	print HTM "<dd>";
	htm_text("$db_subs{$sub_name}{sum}");
    }
    print HTM "</dl>\n";
    print HTM "<HR width=10%>\n";
    print HTM "Depends on the following files: \n";
    print HTM "<dl>\n";
    foreach $f_name ( sort keys %{$db_files{$name}{depend_on_file}} ) {
	print HTM "<dt>";
	htm_link($f_name,"$f_name.html");
	print HTM "<br>";
    }
    print HTM "</dl>\n";
    print HTM "<HR width=10%>\n";
    print HTM "Used by the following files: \n";
    print HTM "<dl>\n";
    foreach $f_name ( @list_db_files ) {
	print HTM "<dt>";
	if( $db_files{$f_name}{depend_on_file}{$name} ) {
	    htm_link($f_name,"$f_name.html");
	    print HTM "<br>";
	}
    }
    print HTM "</dl>\n";
    htm_end();
}


print "printing modules\n";

foreach $name ( keys %db_modules ) {
    htm_start("$output_dir/$name.html","$name");
    print HTM "<table width=\"100%\">\n";
    print HTM "<tr align=center>";
    print HTM "<td><H3>$name</font></H3><td>\n";
    print HTM "<td><H5>File: ";
    htm_link(" $db_modules{$name}{file}","$db_modules{$name}{file}.html");
    print HTM "</H5></td>";
    print HTM "</tr>";
    print HTM "</table>\n";
    #$str = $db_modules{$name}{sum};
    #$str =~ s/>/&gt;/g; $str =~ s/</&lt;/g; 
    #print HTM "Summary: $str<BR>\n";
    htm_text("Summary: $db_modules{$name}{sum}");
    print HTM "Author : $db_modules{$name}{auth}<BR>\n";
    print HTM "Version: $db_modules{$name}{ver}<BR>\n";
    #print HTM "File:";
    #htm_link(" $db_modules{$name}{file}","$db_modules{$name}{file}.html");
    if ( $db_modules{$name}{usage} ) { htm_usage($db_modules{$name}{usage}); }
    if ( $db_modules{$name}{alg} ) { htm_algorithm($db_modules{$name}{alg}); }
    print HTM "<HR width=10%>\n";
    print HTM "Subroutines: \n";
    print HTM "<dl>\n";
    #while (  $sub_name = pop @{$db_modules{$name}{subs}} ) {
    foreach $sub_name ( sort @{$db_modules{$name}{subs}} ) {
	print HTM "<dt>";
	htm_link($sub_name,"$name:$sub_name.html");
	print HTM "<br>";
	$sub_name = "$name:$sub_name";
	#print HTM "<dd>$db_subs{$sub_name}{sum}<BR>\n";
	print HTM "<dd>";
	htm_text("$db_subs{$sub_name}{sum}");
    }
    print HTM "</dl>\n";
    print HTM "<HR width=10%>\n";
    print HTM "Global Variables: \n";
    print HTM "<dl>\n";
    foreach $var ( @{$db_subs{"$name:"}{vars}} ) {
	#$var_name = "$name::$var";
	$var_name = $name.":".":".$var;
	#print "VAR::: $var_name\n";
	if ( $db_vars{$var_name}{decl} =~ /parameter/i ) { $color = "#008800" }
	else { $color = "#880000" }
	if ( ! ($db_vars_used{$var_name} || $db_vars{$var_name}{used_by}) ) {
	    $color = "#ff0000" }
	print HTM "<dt><font color=$color><B>$var</B></font>";
	if ( $db_vars{$var_name}{decl} =~ /^used from (\w+)/ ) {
	    print HTM " : used from "; htm_link("$1",uc("$1").".html");
	}
	else {
	    print HTM " : <code>$db_vars{$var_name}{decl}</code><BR>\n";
	}
	#print HTM "<dd>$db_vars{$var_name}{sum}<BR>\n";
	print HTM "<dd>";
	if ( $db_vars{$var_name}{sum} ) {
	    #print HTM "$db_vars{$var_name}{sum}<BR>\n";
	    htm_text("$db_vars{$var_name}{sum}");
	}
	if ( $db_vars{$var_name}{value} ) {
	  print HTM 
	      "Initial Value <code> = $db_vars{$var_name}{value}</code><BR>";
	}
	if ( ! ($db_vars{$var_name}{sum} || $db_vars{$var_name}{value}) ) {
	    print HTM "<BR>\n";
	}
	#print HTM "      <code>$db_vars{$var_name}{decl}</code><BR>\n";
	my $use_list = $db_vars{$var_name}{used_by};
	if ( $use_list ) {
	    $use_list =~ s/^\s*\|\s*//;
	    #print HTM "Used by: \n";
	    #print HTM "$use_list<BR>\n";
	    print HTM "Used by: \n";
	    my $used_by_sub;
	    foreach $used_by_sub ( sort split / \| /, $use_list ) {
		$used_by_sub =~ s/:$//;
		#print HTM "$used_by_sub<BR>\n";
		print HTM " | ";
		htm_link("$used_by_sub","$used_by_sub.html");
	    }
	    print HTM " |<BR>\n";
	}

    } 
    print HTM "</dl>\n";
    htm_end();
}


print "printing sub's\n";
foreach $name ( keys %db_subs ) {
    htm_start("$output_dir/$name.html","$name");
    #if ( ! $name =~ /([^:]*):([^:]*)/ ) { 
	#print "wrong sub name: $name\n"; next; };
    $sub_name = $name; $sub_name =~ s/^.*://;
    $mod_name = $name; $mod_name =~ s/:.*$//;
    print HTM "<table width=\"100%\">\n";
    print HTM "<tr align=center>";
    print HTM "<td><H3>$sub_name</font></H3><td>\n";
    print HTM "<td><H5>Module: ";
    htm_link("$mod_name","$mod_name.html");
    print HTM "</H5></td>";
    print HTM "<td><H5>File: ";
    htm_link("$db_subs{$name}{file}","$db_subs{$name}{file}.html");
    print HTM "</H5></td>";
    print HTM "</tr>";
    print HTM "</table>\n";


    #print HTM "<H2><Center>$name</Center></font></H2>\n";
    #$str = $db_subs{$name}{sum};
    #$str =~ s/>/&gt;/g; $str =~ s/</&lt;/g; 
    #print HTM "Summary: $str<BR>\n";
    htm_text("Summary: $db_subs{$name}{sum}");
    print HTM "Author : $db_subs{$name}{auth}<BR>\n";
    print HTM "Version: $db_subs{$name}{ver}<BR>\n";
    if ( $db_subs{$name}{usage} ) { htm_usage( $db_subs{$name}{usage} ); }
    if ( $db_subs{$name}{alg} ) { htm_algorithm( $db_subs{$name}{alg} ); }
    print HTM "<HR width=10%>\n";
    print HTM "Declaration:<BR>\n";
    print HTM "<ul><code>$db_subs{$name}{decl}</code></ul>\n";
    if ( $db_subs{$name}{calls} ) {
	print HTM "Calls the following subroutines/functions:<BR>\n";
	print "printing calls for $name\n";
	htm_prt_subs( " ", (split / /, $db_subs{$name}{calls}) );
    }
    print HTM "Variables: \n";
    print HTM "<dl>\n";
    #foreach $var_name ( keys %db_vars ) {
	#if ( $var_name !~ /^$name:(\w+)/ ) { next; }
    #while ( $var = pop @{$db_subs{$name}{vars}} ) {
    foreach $var ( @{$db_subs{$name}{vars}} ) {
	$var_name = "$name:$var";
	#print "VAR::: $var_name\n";
	if ( $db_vars{$var_name}{decl} =~ /parameter/i ) { $color = "#008800" }
	else { $color = "#880000" }
	if ( ! $db_vars_used{$var_name} ) { $color = "#ff0000" }
	print HTM "<dt><font color=$color><B>$var</B></font>";
	if ( $db_vars{$var_name}{decl} =~ /^used from (\w+)/ ) {
	    print HTM " : used from "; htm_link("$1",uc("$1").".html");
	}
	else {
	    print HTM " : <code>$db_vars{$var_name}{decl}</code><BR>\n";
	}
	print HTM "<dd>";
	if ( $db_vars{$var_name}{sum} ) {
	    #print HTM "$db_vars{$var_name}{sum}<BR>\n";
	    htm_text("$db_vars{$var_name}{sum}");
	}
	if ( $db_vars{$var_name}{value} ) {
	    print HTM 
	      "Initial Value <code> = $db_vars{$var_name}{value}</code><BR>\n";
	}
	if ( ! ($db_vars{$var_name}{sum} || $db_vars{$var_name}{value}) ) {
	    print HTM "<BR>\n";
	}
	#print HTM "      <code>$db_vars{$var_name}{decl}</code><BR>\n";
    }
    print HTM "</dl>\n";
    htm_end();

}

print_index();

if ( $component_name =~ /^model$/ ) {
    print_main_index();
}

print "done\n";
print "to view the documentation do: \"netscape $abs_doc_dir/index.html\"\n";


#### this is the end of the main program ######

sub abs_path {
    my $rel_path;
}


sub relpath {
  my ($src,$dst) = @_;
  $src =~ s/[^\/]*\/\.\.\///g;
  $src =~ s/^(?!\/)/\//;
  $src =~ s/\/+$//;
  $dst =~ s/[^\/]*\/\.\.\///g;
  $dst =~ s/^(?!\/)/\//;
  $dst =~ s/\/+$//;
  my $back = 0;
  while(1) {    # relative path
    return '../'x$back.$1 if $dst =~ /$src\/(.*)/;
    $src =~ s/\/[^\/]*$//;
   $back++;
  }
}




sub print_main_index {
    my %run_index = ();
    print  "printing main index\n";
    if ( open ( IND, "$doc_dir/index.html" ) ) {
	while ( <IND> ) {
	    if ( /<!rundeck!><a href=\"([^"]+)\">([^<]+)</ ) { #"))){
	    #if ( /<!rundeck!>/ ) {
		#print "parsing index: $1 $2\n";
		$run_index{$2} = $1;
	    }
	}
    }
    $run_index{ $run_name } = $path_doc_output."/index.html";

    #print "doc, out, rel\n";
    #print "$abs_doc_dir\n";
    #print "$abs_output_dir\n";
    #print relpath($abs_doc_dir,$abs_output_dir)."\n";

    htm_start0("$doc_dir/index.html","GISS GCM Documentation Index");
    print HTM '<H5><Center>GISS GCM Documentation Index</Center></H5>'."\n";
    print HTM "<HR>\n";

    print HTM '<H3><Center>General Documentation</Center></font></H3>'."\n";
#    htm_link("ModelE Reference Manual","modelE.html"); 
#    print HTM "<BR>\n";
#    htm_link("Frequently asked questions about the GISS model","FAQ.html"); 
#    print HTM "<BR>\n";
#    htm_link("HOW-TO document for the GCM","HOWTO.html"); 
#    print HTM "<BR>\n";
#    htm_link("Options for running the GISS GCM","OPTIONS.html"); 
#    print HTM "<BR>\n";

#    while( <$doc_dir/*.txt> ) {
#	s/$doc_dir\///;
#	print "txt loop $_ \n";
#	htm_link("$_", "$_"); print HTM "<BR>\n";
#    }

#    while( <$doc_dir/*.html> ) {
#	s/$doc_dir\///;
#	print "txt loop $_ \n";
#	htm_link("$_", "$_"); print HTM "<BR>\n";
#    }

#    while( <$doc_dir/*/index.html> ) {
    #for 
    foreach $_ ( sort `ls $doc_dir/*/index.html` ) {
	s/$doc_dir\///;
	print "dir loop $_ \n";
	my $doc_link = $_;
	s/\/index.html//;
	htm_link("$_", "$doc_link"); print HTM "<BR>\n";
    }


    print HTM '<H3><Center>Source Code Repository</Center></font></H3>'."\n";
    print HTM 
      "<a href=\"http://simplex.giss.nasa.gov/cgi-bin/gitweb.cgi?p=modelE.git;a=summary\">\n";
    print HTM "View source code in the repository</a>";
    print HTM " for latest updates e.t.c. This link allows you to view \n\
      all the source files currently in Git repository together with their \n\
      older versions. You can also make comparisons between different \n\
      versions of the same file.<BR>\n";
    print HTM "Don't use this link to download the code. Instead read the \n\
      section "; htm_link( " Getting the code from GISS repository",
      "UserGuide/Getting_the_code_form_GISS_repository.html" );
    print HTM " of the "; htm_link( "User Guide",  "UserGuide/index.html");
    print HTM " file.<BR>\n";

    print HTM "<P>\n";
    print HTM '<H3><Center>Dynamic Source Code Documentation</Center></font></H3>'."\n";

    #print HTM "<!rundeck!>";
    #htm_link("new","$output_dir/index.html"); print HTM "<BR>\n";
    
    foreach my $run ( sort keys %run_index ) {
	#print "testing: $run $run_index{$run}\n";
	if ( ! -f "$doc_dir/$run_index{$run}" ) { next; }
	print HTM "Rundeck: \n";
	print HTM "<!rundeck!>";
	htm_link($run,$run_index{$run}); print HTM "<BR>\n";
    }
    htm_end();
 
}

sub print_index {
    print  "printing index\n";
    htm_start("$output_dir/index.html","Rundeck Documentation Index");
    print HTM '<H3><Center>Dynamic Source Code Documentation</Center></font></H3>'."\n";

    htm_link("Source Files","files.html"); print HTM "<BR>\n";
    htm_link("Fortran Modules","modules.html"); print HTM "<BR>\n";
    htm_link("Database Variables","vars_db.html"); print HTM "<BR>\n";
    htm_link("Namelist Variables","vars_nl.html"); print HTM "<BR>\n";
    htm_link("Global Variables ( LONG! )","vars_global.html"); print HTM "<BR>\n";
    htm_link("ALL Variables ( LONG! )","vars_all.html"); print HTM "<BR>\n";

    #print HTM '<H3><Center>Components</Center></font></H3>'."\n";
    print HTM "<BR>\n";
    while( <$output_dir/COMPONENT___*> ) {
	s/$output_dir\///;
	print "found component:  $_ \n";
	htm_link("$_", "$_/index.html"); print HTM "<BR>\n";
    }


    htm_end();
 
}


sub postprocess_lists {

#    print "\n\n>>>>>>>>>before postproc >>>>>>>>>>>>>\n\n\n";
#    foreach $sub_name ( sort keys %db_subs ) {
#	print "$sub_name >>\n";
#	my $use_list = $db_subs{$sub_name}{used_by};
#	print "<< $use_list<<\n\n"; # if ( ! $db_vars{$var_name}{decl} );
#    }

    # make variable declarations look nicer
    foreach $var_name ( keys %db_vars ) {
	$decl = $db_vars{$var_name}{decl};
	$decl =~ s/doubleprecision/real*8/; #make it uniform and shorter
	if ( $decl !~ /$some_decl/i ) { next; } #leave it alone for now
	if ( $decl !~ /^,($some_decl)/i ) {
	    $decl =~ s/,($some_decl)//i;
	    $decl = $1.$decl;
	} else {
	    $decl =~ s/^,//;
	}
	$decl =~ s/,/, /g;
	if ( $decl =~ /^use(\w+)/ ) { 
	    $decl = "used from $1";
	    my $var_from = $1;
	    $var_from =~ tr/a-z/A-Z/;
	    $var_name =~ /^(.*):([^:]*)$/;
	    my $var_prefix = $1; my $var_var= $2;
	    if ( $db_vars{$var_name}{value} =~ />\s*(\w+)/ ) {
		$var_var = $1;  # var renamed with =>
	    }
	    if ( exists $db_vars{"$var_from\:\:$var_var"} ) { #variable
		$db_vars{"$var_from\:\:$var_var"}{used_by} .= " | $var_prefix";
	    } elsif ( exists $db_subs{"$var_from\:$var_var"} ) { #sub
		$db_subs{"$var_from\:$var_var"}{used_by} .= " | $var_prefix";
	    } else { #no object
		print 
		    "USED UNKNOWN: >$var_from\:\:$var_var< in >$var_prefix<\n";
	    }
	}
	$db_vars{$var_name}{decl} = $decl;
    }

#    print "\n\n>>>>>>>>>after postproc >>>>>>>>>>>>>\n\n\n";
#
#    foreach $var_name ( sort keys %db_vars ) {
#	my $use_list = $db_vars{$var_name}{used_by};
#	next if ( ! $use_list );
#	print "$var_name >>\n";
#	print "<< $use_list<<\n\n"; # if ( ! $db_vars{$var_name}{decl} );
#    }
#
#    foreach $sub_name ( sort keys %db_subs ) {
#	my $use_list = $db_subs{$sub_name}{used_by};
#	next if ( ! $use_list );
#	print "$sub_name >>\n";
#	print "<< $use_list<<\n\n"; # if ( ! $db_vars{$var_name}{decl} );
#    }

    # check @calls info and correct it if necessary
    print "postprocessing calls\n";
    foreach $sub_name ( keys %db_subs ) {
	if ( ! $db_subs{$sub_name}{calls} ) { next; } # no call data
	my $calls_str = $db_subs{$sub_name}{calls};
	$calls_str =~ s/\n/ /g;
	$calls_str =~ s/^\s*//;
	$calls_str =~ s/\s*$//;
	#print "CHECKING CALLS: $sub_name\n";
	my @call_list = split /\s*[ ,]\s*/, $calls_str;
	foreach my $call ( @call_list ) {
	    $call =~ tr/A-Z/a-z/;
	    #print "CALLS: $sub_name /  $call\n";
	    if ( $call =~ /:/ ) { 
		($mod_name,$name) = split /:/, $call;
		$call = uc($mod_name).":".lc($name);
		next; 
	    }
	    my $mod_name = (split /:/, $sub_name)[0];
	    #print "MOD_NAME:  $mod_name\n";
	    if ( $db_subs{"$mod_name:$call"}{file} ) {
		$call = "$mod_name:$call";
	    }
	    # just for a test look into global subs:
	    if ( $db_subs{":$call"}{file} ) {
		$call = ":$call";
	    }	    
	}
	$db_subs{$sub_name}{calls} = join " ", @call_list;
	#print ">>> $sub_name >>> $db_subs{$sub_name}{calls}\n";
    }
}

sub by_name_sub_mod { 
    @aa = split /:/,$a; 
    @bb = split /:/,$b; 
    lc($aa[2]) cmp lc($bb[2]) ||  
	lc($aa[1]) cmp lc($bb[1]) || 
	    lc($aa[0]) cmp lc($bb[0])  ;
}

sub by_sub_mod { 
    @aa = split /:/,$a; 
    @bb = split /:/,$b; 
    lc($aa[1]) cmp lc($bb[1]) || 
	lc($aa[0]) cmp lc($bb[0])  ;
}

sub htm_prt_vars { # print list of variables
    my $opts = shift;
    my $var_name;
    my $file;
    my $print_location = $opts =~ /l/i;
    my $do_sort        = $opts =~ /s/i;
    my @var_list;
    if ( $do_sort ) {
	@var_list = sort by_name_sub_mod @_;
    } else {
	@var_list = @_;
    }
    my $first_letter;
    foreach $first_letter ( 'a' .. 'z' ) {
	print HTM "  <A HREF=\"#$first_letter\">$first_letter</A>";
    }
    print HTM "<dl>\n";
    $first_letter = "a";
    print HTM "<A NAME=a></A>\n";
    foreach $var_name ( @var_list ) {
	$var_name =~ /(\w*):(\w*):(\w+)/ || next;
	my $mod = $1;
	my $subr = $2;
	my $var = $3;
	#print "PRT_VARS: $var_name, $mod, $subr, $var\n";
	if ( $var !~ /^$first_letter/ ) {
	    my $first_letter_prev = ++$first_letter;
	    $var =~ /^(.)/; $first_letter = $1;
	    foreach my $letter ( $first_letter_prev .. $first_letter ) {
		print HTM "<A NAME=$letter></A>\n";
	    }
	}
	if ( $db_vars{$var_name}{decl} =~ /parameter/i ) { $color = "#008800" }
	else { $color = "#880000" }
	if ( ! $db_vars_used{$var_name} ) { $color = "#ff0000" }
	print HTM "<dt><font color=$color><B>$var</B></font>";
	print HTM " : <code>$db_vars{$var_name}{decl}</code><BR>\n";
	print HTM "<dd>";
	if ( $print_location ) {
	    if ( $subr ) {
		#print HTM "Subroutine: $subr, ";
		print HTM "Subroutine: ";
		htm_link("$subr","$mod:$subr.html");
		$file = $db_subs{"$mod:$subr"}{file};
	    } else {
		print HTM "Global variable ";
		$file = $db_modules{"$mod"}{file};
	    }
	    #print HTM "Module: $mod, File: $file<BR>\n";
	    print HTM ". Module: ";
	    if ( $mod ) { 
		htm_link("$mod","$mod.html");
	    } else {
		print HTM "NONE";
	    }
	    print HTM ". File: ";
	    htm_link("$file","$file.html");
	    print HTM "<BR>\n";
	}
	if ( $db_vars{$var_name}{sum} ) {
	    htm_text("$db_vars{$var_name}{sum}");
	}
	if ( $db_vars{$var_name}{value} ) {
	    print HTM 
	      "Initial Value <code> = $db_vars{$var_name}{value}</code><BR>\n";
	}
	if ( ! ($db_vars{$var_name}{sum} || $db_vars{$var_name}{value}) ) {
	    print HTM "<BR>\n";
	}
    }
    if ( $first_letter !~ /z/ ) {
	++$first_letter;
	foreach my $letter ( $first_letter .. 'z' ) {
	    print HTM "<A NAME=$letter></A>\n";
	}
    }
    print HTM "</dl>\n";

}


sub htm_prt_subs {
    my $opts = shift;
    my $sub_name;
    my $file;
    my $print_location = $opts =~ /l/i;
    my $do_sort        = $opts =~ /s/i;
    my @sub_list;
    if ( $do_sort ) {
	@sub_list = sort by_sub_mod @_;
    } else {
	@sub_list = @_;
    }

    #print "printing SUBS: ", join " ", @sub_list, "\n"; 

    print HTM "<dl>\n";
    foreach $sub_name ( @sub_list ) {
	if ( $db_subs{$sub_name}{file} ) {
	    my $file = $db_subs{$sub_name}{file};
	    my $name, $mod_name;
	    ($mod_name,$name) = split /:/, $sub_name;
	    print HTM "<dt>";
	    htm_link($name,"$sub_name.html");
	    print HTM "<dd>";
	    print HTM "Module: "; 
	    if ( $mod_name ) { htm_link($mod_name,"$mod_name.html"); }
	    else { print HTM "Global "; }
	    print HTM " File: "; htm_link($file,"$file.html");
	    print HTM "<BR>\n";
	    htm_text("$db_subs{$sub_name}{sum}");
	} else {
	    print "Warning: Can't find \"$sub_name\" listed in \@calls\n";
	    print "    Be more specific, use \@calls module_name:sub_name\n";
	    print HTM "<dt>$sub_name";
	    print HTM "<dd><br>\n";
	}
    }
    print HTM "</dl>\n";
}

sub htm_start0 {
    my $filename = shift;
    my $title = shift;
    open ( HTM, ">$filename" ) || die "can't open $filename";
    print HTM "<html>\n";
    print HTM "<head>\n";
    print HTM "<TITLE>$title</TITLE>\n";
    print HTM "</head>\n";
    print HTM '<body BGCOLOR="#FFFFFF"TEXT="#000000">'."\n";
    if ($html_include_dir) {
	print HTM "<!--#include virtual=\"$html_include_dir/modelE_header.html\" -->\n";
    }
}

sub htm_start {
    my $filename = shift;
    my $title = shift;
    htm_start0( $filename, $title );
    print HTM "<table width=\"100%\">\n";
    print HTM "<tr align=left>";
    print HTM "<td>";
    htm_link("Index","$path_output_doc/index.html");
    print HTM "</td><td>";
    print HTM "Rundeck: ";
    htm_link($run_name,"index.html");
    print HTM "</td><td align=right>";
    print HTM "<i>Created: $current_date</i>";
    print HTM "</td>";
    print HTM "</tr>";
    print HTM "</table>\n";
}

sub htm_end {
    if ($html_include_dir) {
	print HTM "<!--#include virtual=\"$html_include_dir/footer.html\" -->\n";
    }
    print HTM "</body>\n";
    print HTM "</html>\n";
    close HTM;
}

sub htm_link {
    my $name = shift;
    my $link = shift;
    $link =~ s/:/%3A/g;
    print HTM "<a href=\"$link\">$name</a>";
}

sub htm_text {
    my $str = shift;
    $str =~ s/&/&amp;/g; $str =~ s/>/&gt;/g; $str =~ s/</&lt;/g; 
    print HTM "$str<BR>\n";
}

sub htm_usage {
    my $str = shift;
    print HTM "Usage:<BR>\n";
    print HTM "<ul><pre>\n";
    htm_text($str);
    print HTM "</pre></ul>\n";
}

sub htm_algorithm {
    my $str = shift;
    print HTM "Algorithm:<BR>\n";
    print HTM "<ul><pre>\n";
    htm_text($str);
    print HTM "</pre></ul>\n";
}

sub htm_incl_file {
    print HTM "<pre>\n";
    while ( <DOCFILE> ) {
	s/&/&amp;/g; s/>/&gt;/g; s/</&lt;/g;
	print HTM;
    }
    print HTM "</pre>\n";
}

#@sum   UNITNAME Brief summary/description of the current program unit or file
#@auth    Author  
#@ver     Version (a number plus any relevant comments)
#@calls   List of routines called in this program unit
#@cont    List of program units contained within a file or module + other info.
#@fun     FUNCNAME Denotes a function
#@param   PARAMNAME Denotes a parameter (in the FORTRAN sense)
#@var     VARNAME Denotes a variable (in the FORTRAN sense)
#@dbparam VARNAME Denotes a database parameter (in the modelE sense)
#@nlparam VARNAME Denotes a NAMELIST parameter (in the modelE sense)
#@+       Continuation line for @sum/@calls/@cont

sub parse_file {
    my $free_form = 0;
    if ( $current_file =~ /\.F90$/ ) { $free_form = 1; }
    #print "parsing $current_file\n";
    # resetting globals
    $current_module = "";
    $current_sub = "";
    $current_typedef = "";
    while( <SRCFILE> ) {
	chop;
	#strip regular comments
	if ( (/^C/i && !$free_form) || /^![^@]/ ) { next; }
	#parse !@... info here
#	if( /^!\@sum/i ) { #subroutine summary
#	    $doc_tag = "sum";
#	    s/^!\@sum\s+(\w+)\s*//i;
#	    $sum_name = "$current_module:$1";
#	    $current_sub = $1;
#	    $db_subs{$sum_name}{sum} = "$_\n";
#	    next;
#	}

	#strip spaces at the end
	s/\s*$//;

	if( /^!\@(var|param|dbparam|nlparam)\b/i ) { #variable summary
	    $doc_tag = "var";
	    s/^!\@(\w+)\s+((\w+)(\s*,\s*(\w+))*)\s*//i;
	    $tag = $1;
	    #@var_list = split '[ ,]', $1;
	    @var_list = split /\s*,\s*/, $2;
	    while( $name = pop @var_list ) {
		$name =~ tr/A-Z/a-z/;
		$var_name = "$current_module:$current_sub:$name";
		#print "var: $var_name\n";
		#print "comment: $_\n";
		$db_vars{$var_name}{sum} = "$_\n";
		$db_vars{$var_name}{tag} = $tag;
		push @{$db_subs{"$current_module:$current_sub"}{vars}}, $name;
	    }
	    next;
	}

	if( /^!\@(\w+)\b/i ) { #generic tag
	    $doc_tag = $1;
	    s/^!\@\w+\s*//i;
	    if ( $current_sub ) {
		$db_subs{"$current_module:$current_sub"}{$doc_tag} = "$_\n";
		    #print "Sub tag:\n";
		    #print "$current_sub  $doc_tag : $_\n";
	    } elsif ( $current_module ) {
		$db_modules{"$current_module"}{$doc_tag} = "$_\n";
		    #print "Module tag:\n";
		    #print "$current_module  $doc_tag : $_\n";
	    } else {
		$db_files{"$current_file"}{$doc_tag} = "$_\n";
		    #print "File tag:\n";
		    #print "$current_file  $doc_tag : $_\n";
	    }
	    next;
	}

	if( /^!\@\+/i ) { #continuation line
	    #$doc_tag = "var";
	    s/^!\@\+\s?\s?\s?//i;
	    if ( $doc_tag =~ "var" ) {
		while( $name = pop @var_list ) {
		    $var_name = "$current_module:$current_sub:$name";
		    $db_vars{$var_name}{sum} .= "$_\n"; }
	    #} elsif ( $doc_tag =~ "sum" ) {
		#$db_subs{$sun_name}{sum} .= "$_\n";
	    } else { #generic
		if ( $current_sub ) {
		    $db_subs{"$current_module:$current_sub"}{$doc_tag}.="$_\n";
		} elsif ( $current_module ) {
		    $db_modules{"$current_module"}{$doc_tag} .= "$_\n";
		    #print "Module tag:\n";
		    #print "$current_module  $doc_tag : $_\n";
		} else {
		    $db_files{"$current_file"}{$doc_tag} .= "$_\n";
		}
	    }
	    next;
	}
	#print;
	if ( /hjdhsdjhfsjkdf/ ) { print "you are kidding!\n"; }

	if( /^.+!\@(var|param|dbparam|nlparam)\b/i ) { #inline variable summary
	    $doc_tag = "var";
	    $str = $_;
	    $str =~ s/.*!\@(\w+)\s+((\w+)(\s*,\s*(\w+))*)\s*//i;
	    $tag = $1;
	    #$str =~ s/\s*!\@var\s+(.*)$//;
	    #$comment = $1;
	    #$str =~ s/......//;
	    #$str =~ /(
		#      (\w+ (\s*=\s*[^,]+)? \s*,\s*)* \w+ (\s*=\s*[^,]+)?
		#      )\s*$/x;
	    @var_list = split /\s*,\s*/, $2;
	    while( $name = pop @var_list ) {
		$name =~ tr/A-Z/a-z/;
		$var_name = "$current_module:$current_sub:$name";
		$db_vars{$var_name}{sum} = "$str\n";
		$db_vars{$var_name}{tag} = $tag;
		push @{$db_subs{"$current_module:$current_sub"}{vars}}, $name;
	    }
	    s/!\@.*$//i;
	}

	# prarse fortran code here
	if ( /!/ ) {  # possible comment at the en of the line
	    /^(
	       (  [^"'!] | (\"[^"]*\") | (\'[^']*\')  )*
              )/x;  # '))) / ; -emacs fix
	    $_ = $1;
	    #print "$_\n";
	}

	if ( $free_form ) {
	    if ( $fstr =~ /\\$/ ) { # continuation line
		 chop $fstr;
		 s/^\s*\\//;
		 $fstr .= $_;
		next;
		}
	} else {
	    if (  /^     \S/ || /^\s*$/ ) { # continuation line
		s/^     \S/ /;
		$fstr .= $_;
		next;
	    }
	}

	# parse previous line
	parse_fort_str();

	#first check if beginning/end of block

	if ( /^\s*module\s+(\w+)/i ) { #start module
	    if ( $1 !~ /procedure/i ) {
		$current_module = uc($1);
		$current_sub = "";
		$db_modules{"$current_module"}{file} = $current_file;
		push @{$db_files{"$current_file"}{modules}}, $current_module;
		#print "MOD: ", ${$db_files{"$current_file"}{modules}}[0],"\n";
	    }
	}
	if ( /^\s*end\s+module\b/i ) { #end module
	    $current_module = "";
	}
	if ( 
/^\s*(subroutine|(?:$some_decl)?\s*function|program|interface)\s+(\w+)/i ) {
	    $current_sub = lc($2);
	    $db_subs{"$current_module:$current_sub"}{file} = $current_file;
	    if ( $current_module ) {
		push @{$db_modules{"$current_module"}{subs}}, $current_sub;
	    } else {
		push @{$db_files{"$current_file"}{subs}}, $current_sub;
		#push @{$db_modules{"\@GLOBAL"}{subs}}, $current_sub;
	    }
	}
	if ( /^\s*end\s+(subroutine|function|program|interface)\b/i 
	     || /^\s*end\s*$/i ) {
	    $current_sub = "";
	}
	if ( /^\s*(block\s+data)\s+(\w+)/i ) { # block data hack
	    $current_sub = 'BLOCK_DATA_'.$2;
	}
	if ( /^\s*type\s+(\w+)/i ) {
	    $current_typedef = uc($1);
	}
	if ( /^\s*end\s+type\b/i ) {
	    $current_typedef = "";
	}

	# beginning of the new fortran line
	$fstr = $_;
				   
    }
    print "number of lines: $.\n";
}

sub parse_fort_str {
    #print "$fstr\n";

    # subroutine/function declaration
    if ( $fstr =~ /^\s*(subroutine|(?:$some_decl)?\s*function)/i ) {
	my $decl = lc($fstr);
	$decl =~ s/^\s*//; $decl =~ s/\s*$//;
	$decl =~ s/\s*,\s*/, /g;
	$db_subs{"$current_module:$current_sub"}{decl} = $decl;
	return;
    }

    if ( $current_typedef ) { return; } #skip typedefs for now

#!!! experimental code
# create the list of USEd modules
    if ( $fstr =~ /^\s*use\s+(\w+)/i ) {
	$db_files{"$current_file"}{depend_on_mod}{uc($1)} = 1;
    }
  
    # variable declaration string
    if ( $fstr =~ s/^\s*($some_decl)\s*//i ) {
	$var_type = $1;
	if ( $fstr =~ s/([^"']*)\s*::\s*//i ) {        #"))){
	    $var_type .= $1;
	}
    #if ( $fstr =~ /^\s*(integer|real|double|character|logical|complex)/i ) {
	#if ( $fstr =~ s/^\s*(\S+.*)::\s*//i ) {
	#    $var_type = $1;
	#} else {
	#    $fstr =~ s/^\s*(\S+)\s*//i;
	#    $var_type = $1;
	#    if ( $var_type =~ /double/i ) {
	#	$var_type = "real*8";
	#	$fstr =~ s/^\s*(\S+)\s*//i;
	#    }
	#}
	$var_type =~ s/ //g;
	$var_type =~ tr/A-Z/a-z/;
	#print "VAR_TYPE:  $var_type\n";
	# parse single variables
	while ( $fstr =~ s/^\s*,*\s*(\w+)((\([^()]+\))*)// ) {
	    $var = $1;
	    $dim = $2;
	    $var =~ tr/A-Z/a-z/;
	    #print "VAR $var   DIM $dim\n";
	    # look for initial value
	    $var_val = "";
	    # var = abc( ()), ( ())
	    if( $fstr =~ 
		s/^\s*=\s*(
			   [^,"'()]* (   \(  ( [^()] | (\([^()]*\)) )*  \)  )
                          )//x ) {
              #  " ))) {
		$var_val = " $1";
		#print "VAR:  $var   ::: $var_val\n";
	    }
	    # var = "abc", 'abc', 123
	    if( $fstr =~ s/^\s*=\s*(\"[^"]*\"|\'[^']*\'|[^,"'()]+)// ) {#'))) {
		#print "FSTR: #$fstr#\n";
		#if( $fstr =~ s/^\s*=\s*(\".+\")// ) {
		#print "VAL: #$1#\n";
		$var_val = " $1";
		#print "VAR:  $cvar   ::: $var_val\n";
	    }
	    $var_name = "$current_module:$current_sub:$var";
	    $db_vars_declared{$var_name} = 1;
	    $db_vars{$var_name}{decl} .= ",$var_type";
	    if ( $dim ) { $db_vars{$var_name}{decl} .= ",dimension$dim"; }
	    $db_vars{$var_name}{value} = $var_val;
	    if ( ! $db_vars{$var_name}{sum} ) { #wasn't included earlier
		push @{$db_subs{"$current_module:$current_sub"}{vars}}, $var;
	    }
	    #if ( ! $current_sub ) {
		#print "GLOBAL: $current_module $var_name\n";
	    #}
	}
	return;
    }
    
    # parse some extra declarations
    if ( $fstr =~ s/^\s*dimension\s*//i ) {
	while ( $fstr =~ s/^\s*,*\s*(\w+)\s*($parenth3)//i ) {
	    $var = $1;
	    $dim = $2;
	    $var =~ tr/A-Z/a-z/;
	    $var_name = "$current_module:$current_sub:$var";
	    $db_vars{$var_name}{decl} .= ",dimension$dim";
	}
	return;
    }

    # skip common blocks
    if ( $fstr =~ /^\s*(common)/i ) { 
	return;
    }

#!!! experimental code
# check for unused variables (very primitive check)
    $tmp_fstr = $fstr;
    $tmp_fstr =~ s/\".*?\"//g; $tmp_fstr =~ s/\'.*?\'//g;
    foreach $var ( split /\W+/, $tmp_fstr ) {
	$var =~ tr/A-Z/a-z/;
	$var_name = "$current_module:$current_sub:$var";
	if ( $db_vars_declared{$var_name} ) {
	    $db_vars_used{$var_name} = 1;
	} else {
	    $var_name = "$current_module\:\:$var";
	    if ( $db_vars_declared{$var_name} ) {
		$db_vars_used{$var_name} = 1;
	    }
	}
    }
}

sub old_old {
while ($filename = shift) {
  open(IN, $filename) or die "can't open $filename\n";
  #print "reading $filename\n";
  $str = "";
  $function = "_GLOBAL_";
  $module = "_NONE_";
  $print_header = 1;
  while (<IN>) {
    if ( /^[Cc]/ || /^ *!/ || /^\s*$/ ) {next;}
    if (  /^     \S/ ) {
      #print "cont: $_\n";
      $str .= $_;
      #print;
    }
    else {
      if ( $str =~ /^\s*(module)\s+(\w+)/i ) {
        $module = $2;
        $function = "_GLOBAL_";
        #print "FUCTTION:: $function\n";
        $print_header = 1;
      }
      if ( $str =~ /^\s*(subroutine|function|program)\s+(\w+)/i ) {
        $function = $2;
        #print "FUCTTION:: $function\n";
        $print_header = 1;
      }
      if ( $str =~ /^\s*(block\s+data)\s+(\w+)/i ) {
        $function = 'BLOCK_DATA_'.$2;
        #print "FUCTTION:: $function\n";
        $print_header = 1;
      }
      if ( $str =~ /\b$var\b/i ) { 
        $is_decl = ($str =~ /^\s*(integer|real|double|logical|character|data|dimension)/i);
        $is_used = ($str =~ /^\s*(use|common)/i);
        $is_changed = ($str =~ /\b$var\b\s*(\(.*\))?\s*=|call.*\b$var\b/i);
        if( ($show_decl && $is_decl) || 
            ($show_use && $is_used) ||
            ($show_body && ( (!($is_decl||$is_used)) || $str=~/=> *\b$var\b/i))
          ||($show_changed && $is_changed) ) {
          if ( $print_header ) {
            print "\nFILE:: $filename  MODULE:: $module  SUB:: $function\n\n";
            $print_header = 0;
          }
          $str =~ s/\b($var)\b/&to_bold($1)/egi;
          print "$. :\n";
          print $str;
        }
      }
      $str = $_;
    }
  }
}

} # sub old_old

sub to_bold {
  my $str = shift(@_);
  #print "got $str\n";
  $str =~ s/(.)/$1\x08$1/g;
  #print "made: $str\n";
  return $str;
}

