#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    script=`basename @script@`
    cat $script.log | sed -e 's/^| Time estimate.*$/| Time estimate [...]/' > $script.log.tmp && mv $script.log.tmp $script.log

    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
