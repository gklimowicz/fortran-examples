#!/usr/bin/perl

# To skip dependencies on certan modules list them here
# for example:
# $skip_modules="domain_decomp.mod|model_com.mod|constant.mod";
# In particuler recent version of the Fortran standard (2003 and 2008)
# have introduced the concept of "intrinsic modules" that are provided
# by the compiler.
$skip_modules="iso_c_binding.mod|iso_fortran_env.mod|ieee_arithmetic.mod|pfunit_mod.mod";

while(<>) {

    # hack to skip unneeded dependencies
    if ( $skip_modules ) { s/ ($skip_modules)//g; }

    $MOD = "mod"; # could be different for other compilers

    if ( /\s*(\w+)\.$MOD:\s+(\w+)\.o\s*$/ ) {
	print "$1\@$2.smod: $2.o\n";
	print "$1.$MOD: $1\@$2.smod\n";
    } else {
	print ;
    }


}
