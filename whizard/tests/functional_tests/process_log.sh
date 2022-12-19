#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --no-model
    cat ${s}_1_p1.log | sed -e 's/^   T(10k evt) =.*$/   T(10k evt) =  <time estimate>/'  -e 's/^\$fc =.*$/\$fc = Fortran-compiler/' -e 's/^\$fcflags =.*$/\$fcflags = Fortran-flags/' -e 's/^\$fclibs =.*$/\$fclibs = Fortran-libs/' -e 's/^\?omega_openmp =.*$/\?omega_openmp => false/' -e 's/^\?openmp_is_active =.*$/\?openmp_is_active = false/' -e 's/^openmp_num_threads_default =.*$/openmp_num_threads_default = 1/' -e 's/^real_range =.*$/real_range = <real_range>/' -e 's/^real_precision =.*$/real_precision = <real_precision>/' -e 's/^real_epsilon =.*$/real_epsilon = <real_epsilon>/' -e 's/^real_tiny =.*$/real_tiny = <real_tiny>/' > ${s}_1_p1.log.tmp
    mv ${s}_1_p1.log.tmp ${s}_1_p1.log    
    diff ref-output/`basename @script@`.ref `basename @script@`_1_p1.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
