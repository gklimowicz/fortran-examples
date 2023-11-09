#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
script=`basename @script@`
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    echo "Contents of ${script}p.evt:" >> $script.log
    cat ${script}p.evt >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
