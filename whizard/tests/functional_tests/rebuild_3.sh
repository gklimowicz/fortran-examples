#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --no-model
    mv $script.log $script.1.log
    ./run_whizard.sh @script@ --no-logging --no-model
    mv $script.log $script.2.log
    cat $script.1.log $script.2.log > $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
