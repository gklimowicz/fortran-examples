#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging --no-model
    script=`basename @script@`
    echo "Unmodified BR:" >> $script.log
    grep "branching ratio" ${script}_zh_a.debug >> $script.log
    echo "Modified BR:" >> $script.log
    grep "branching ratio" ${script}_zh_b.debug >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
