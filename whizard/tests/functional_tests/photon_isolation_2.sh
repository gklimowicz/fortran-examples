#!/bin/sh
### Check WHIZARD for photon isolation from EM acitivity
echo "Running script $0"
script=`basename @script@`
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
