#!/bin/sh
# Testing separate NLO calculation for the Born,
# the virtual and the real component of ee -> jj
# and the simulation of events for each component
echo "Running script $0"
if test -f OCAML_FLAG -a -f FASTJET_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    echo "Contents of ${name}_p1.debug:" >> $name.log
    cat ${name}_p1.debug >> $name.log
    echo "Contents of ${name}_p2.debug:" >> $name.log
    cat ${name}_p2.debug >> $name.log
    cat ${name}_p2_fks_regions.out >> $name.log
    echo "Contents of ${name}_p3.debug:" >> $name.log
    cat ${name}_p3.debug >> $name.log
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements / FastJet available, test skipped"
    exit 77
fi
