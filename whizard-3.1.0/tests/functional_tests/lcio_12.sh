#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f LCIO_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No LCIO or O'Mega matrix elements available, test skipped"
    exit 77
fi
