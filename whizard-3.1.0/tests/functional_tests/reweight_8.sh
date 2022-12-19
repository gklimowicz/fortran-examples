#!/bin/sh
### Check WHIZARD for a simple test process
name=`basename @script@`
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging --no-model
    echo "Contents of ${name}a.weights.dat" >> $name.log
    cat ${name}a.weights.dat >> $name.log
    echo "Contents of ${name}b.weights.dat" >> $name.log
    cat ${name}b.weights.dat >> $name.log
    echo "Contents of ${name}c.weights.dat" >> $name.log
    cat ${name}c.weights.dat >> $name.log
    echo "Contents of ${name}d.weights.dat" >> $name.log
    cat ${name}d.weights.dat >> $name.log
    diff ref-output/`basename @script@`.ref `basename @script@`.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
