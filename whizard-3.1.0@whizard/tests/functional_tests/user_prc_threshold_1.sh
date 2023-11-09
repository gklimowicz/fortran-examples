#!/bin/sh
### Check WHIZARD: process library handling
echo "Running script $0"
if test -f OCAML_FLAG; then
    name=`basename @script@`
    rm -f SM_tt_threshold.grid
    ./run_whizard.sh @script@ --no-logging
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
