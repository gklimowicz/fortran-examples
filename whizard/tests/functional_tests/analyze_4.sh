#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f FASTJET_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    script=`basename @script@`
    echo "Contents of file $script.dat:" >> $script.log
    cat $script.dat >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements / FastJet available, test skipped"
    exit 77
fi
    
