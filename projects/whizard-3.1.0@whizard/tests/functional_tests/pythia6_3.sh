#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f PYTHIA6_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    script=`basename @script@`
    echo "Partial contents of ${script}a.debug:" >> $script.log
    cat ${script}a.debug >> $script.log
    echo "Partial contents of ${script}b.debug:" >> $script.log
    cat ${script}b.debug >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available and/or PYTHIA6 disabled, test skipped"
    exit 77
fi
