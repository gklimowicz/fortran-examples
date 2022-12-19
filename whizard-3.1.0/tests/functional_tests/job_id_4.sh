#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    script=`basename @script@`
    rm -rf ${script}_x*
    ./run_whizard.sh @script@ --no-logging --no-model -J "8001"
    echo "* Files created by integrate:" >> $script.log
    ls ${script}_x.*.1.* >> $script.log
    echo "* Files created in integrate workspace:" >> $script.log
    ls ${script}_x.*/* >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
