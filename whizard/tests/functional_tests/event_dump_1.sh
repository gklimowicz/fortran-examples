#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
script=`basename @script@`
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging --no-model
    echo "Contents of ${script}a.pset.dat:" >> $script.log
    cat ${script}a.pset.dat >> $script.log
    echo "Contents of ${script}b.pset.dat:" >> $script.log
    cat ${script}b.pset.dat >> $script.log
    echo "Contents of ${script}c.pset.dat:" >> $script.log
    cat ${script}c.pset.dat >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
