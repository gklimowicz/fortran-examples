#!/bin/bash
# eda.sh -- Use m4 and a Doxyfile to expand the Doxygen aliases in the documentation.
filter() { sed 's/\\rsb/RSB_D_/g;s/\\librsb\>/RSB_D_librsb/g;s/\\see_/RSB_D_see/g'| sed "s/#/"'\`'"#'/g"; }
cat doc/Doxyfile | grep -v "'" | grep -v '^#' | grep ^ALIASES | sed 's/^[^"]\+//g' | sed 's/^"/\\/g;s/" *$//g' |  filter | sed 's/\([^=]\+\)=\([^=]*\)/define(\`\1'"'"', \`\2'"'"')dnl/g'   > Doxyfile.m4
( echo "include("'`'"Doxyfile.m4')dnl " ; cat rsb_rsb.c ; ) | filter > rsb_rsb.m4
cat rsb_rsb.m4 | m4 -I . > rsb_rsb_e.c
grep RSB_D_ rsb_rsb_e.c && false
