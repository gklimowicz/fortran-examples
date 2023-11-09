#!/bin/sh
### Check WHIZARD for a simple test process
### (This test has no ref file, it just compares event output data)
echo "Running script $0"
if test -f OCAML_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model QCD
    echo "Contents of ${name}.lhe" >> $name.log
    # Maybe we want to compare to reference output data though the shower is
    # quite sensitive...
    #cat ${name}.lhe >> $name.log
    #diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
