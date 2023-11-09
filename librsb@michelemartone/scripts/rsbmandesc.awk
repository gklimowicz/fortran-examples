#!/usr/bin/awk -f
#
# This program prints out the Grammar section of a yacc.output file. 
#
BEGIN { sp=0; }
/^.SH NAME$/ { sp=sp+1; }
/.*/ {
		if(sp==1){sp=sp+1;}
		else
		if(sp==2)
		{
			print "librsb - ";
			sp=sp+1;
		}
		else
		if(sp==3)
		{
			print ".SH DESCRIPTION";
			sp=sp+1;
		}
		print;
	}
