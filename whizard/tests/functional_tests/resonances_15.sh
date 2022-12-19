#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f LCIO_FLAG -a -f PYTHIA6_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    echo "Checking consistency of event weights:" >> $script.log
    ./${script}_check >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No LCIO or O'Mega matrix elements or Pythia6 available, test skipped"
    exit 77
fi
