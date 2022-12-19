#!/bin/sh
### Check WHIZARD for a simple test process
name=`basename @script@`
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging --no-model
    echo "Contents of ${name}_p1.weights.dat" >> $name.log
    cat ${name}_p1.weights.dat >> $name.log
    diff ref-output/`basename @script@`.ref `basename @script@`.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
