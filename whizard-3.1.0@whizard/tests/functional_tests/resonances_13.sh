#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    echo "Output of ${script}_a.pset.dat:" >> $script.log
    cat ${script}_a.pset.dat >> $script.log
    echo "Output of ${script}_b.pset.dat:" >> $script.log
    cat ${script}_b.pset.dat >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
