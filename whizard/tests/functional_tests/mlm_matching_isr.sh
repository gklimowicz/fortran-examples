#!/bin/sh
### Check WHIZARDs initial state shower/matching
echo "Running script $0"
if test -f OCAML_FLAG -a -f PYTHIA6_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    script=`basename @script@`
    echo "Contents of ${script}_a.debug:" >> $script.log
    grep momenta ${script}_a.debug >> $script.log
    echo "Contents of ${script}_b.debug:" >> $script.log
    grep momenta ${script}_b.debug >> $script.log
    diff ref-output/$script.ref $script.log

else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements and/or PYTHIA6 disabled, test skipped"
    exit 77
fi
