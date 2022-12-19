#!/bin/sh
### Check WHIZARD with its variable setup:
echo "Running script $0"
s=`basename @script@`
./run_whizard.sh @script@ --no-logging
cat $s.log | sed -e 's/^\$fc =>.*$/\$fc => Fortran-compiler/' -e 's/^\$fcflags =>.*$/\$fcflags => Fortran-flags/' -e 's/^\$fclibs =>.*$/\$fclibs => Fortran-libs/' -e 's/^\?omega_openmp =.*$/\?omega_openmp => false/' -e 's/^\?openmp_is_active\* =.*$/\?openmp_is_active\* = false/' -e 's/^openmp_num_threads_default\* =.*$/openmp_num_threads_default\* = 1/' > $s.log.tmp
mv $s.log.tmp $s.log
diff ref-output/`basename @script@`.ref `basename @script@`.log
