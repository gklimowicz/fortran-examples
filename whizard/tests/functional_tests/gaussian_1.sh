#!/bin/sh
### Check WHIZARD/O'Mega gaussian-beam spread setup
echo "Running script $0"
if test -f OCAML_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model SM
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
