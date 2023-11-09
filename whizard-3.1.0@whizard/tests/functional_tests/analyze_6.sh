#!/bin/sh
### Check WHIZARD for a simple test process
name=`basename @script@`
echo "Running script $0"
if test -f OCAML_FLAG; then
    rm -f ${name}.dat
    ./run_whizard.sh @script@ --no-logging --no-model
    ./${name}_check
    diff ref-output/$name.ref ${name}_check.out
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
