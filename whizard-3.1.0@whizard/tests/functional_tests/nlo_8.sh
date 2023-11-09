#!/bin/sh
# Testing combined NLO calculation of ee -> jj
# as well as the simulation of combined events
echo "Running script $0"
if test -f OCAML_FLAG -a -f FASTJET_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    echo "Contents of ${name}_p1.debug:" >> $name.log
    cat ${name}_p1_fks_regions.out >> $name.log
    cat ${name}_p1.debug >> $name.log
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements / FastJet available, test skipped"
    exit 77
fi
