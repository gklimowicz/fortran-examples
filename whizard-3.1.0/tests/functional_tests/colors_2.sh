#!/bin/sh
### Check WHIZARD/O'Mega with color correlations
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    script=`basename @script@`
    echo "Events for fixed-energy beams:" >> $script.log
    cat ${script}_fix.pset.dat >> $script.log
    echo "Events for beams with PDF:" >> $script.log
    cat ${script}_pdf.pset.dat >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
