#!/bin/sh
### Check WHIZARD/O'Mega with externally linked LHAPDF library
echo "Running script $0"
if test -f OCAML_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model SM
    cat ${script}a.pset.dat >> $script.log
    cat ${script}.obs.dat >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
