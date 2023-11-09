#!/bin/sh
### Check WHIZARD rerunning with model change
echo "Running script $0"
if test -f OCAML_FLAG; then
    rm -f @script@_lib.* @script@_p?.*
    s=`basename @script@`
    ./run_whizard.sh @script@ --model SM --no-logging
    cat ${s}_p_i1.f90 > $s.1.f90
    ./run_whizard.sh @script@ --no-rebuild --model QED --no-logging
    cat ${s}_p_i1.f90 > $s.2.f90
    # Files should be different
    ! diff -sq $s.1.f90 $s.2.f90
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
