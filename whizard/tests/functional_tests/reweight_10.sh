#!/bin/sh
### Check WHIZARD for a simple test process
name=`basename @script@`
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging --no-model
    cat $name.log | sed -e 's/WHIZARD 3.*$/WHIZARD [version]/' > $name.log.tmp
    mv $name.log.tmp $name.log
    echo "Contents of ${name}_in.weights.dat" >> $name.log
    cat ${name}_in.weights.dat >> $name.log
    echo "Contents of ${name}_out.weights.dat" >> $name.log
    cat ${name}_out.weights.dat >> $name.log
    diff ref-output/`basename @script@`.ref `basename @script@`.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
