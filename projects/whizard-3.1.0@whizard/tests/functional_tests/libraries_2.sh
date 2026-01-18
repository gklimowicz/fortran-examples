#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    rm -f @script@_lib.* @script@_p?.*
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    cat $s.log > $s.1.log
    ./run_whizard.sh @script@ --no-rebuild --no-logging
    cat $s.log > $s.2.log
    ./run_whizard.sh @script@ --model QED --no-rebuild --no-logging
    cat $s.log > $s.3.log
    cat $s.1.log $s.2.log $s.3.log > $s.log
    diff ref-output/$s.ref $s.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    

