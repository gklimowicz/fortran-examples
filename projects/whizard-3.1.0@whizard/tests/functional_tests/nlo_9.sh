#!/bin/sh
# Testing separate NLO calculation of pp -> ee
# as well as the simulation of events for all components
# using dummy-output for the matrix elements
echo "Running script $0"
if test -f OCAML_FLAG ; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    echo "Contents of ${name}_p1.debug:" >> $name.log
    cat ${name}_p1.debug >> $name.log
    echo "Contents of ${name}_p2.debug:" >> $name.log
    cat ${name}_p2_fks_regions.out >> $name.log
    cat ${name}_p2.debug >> $name.log
    echo "Contents of ${name}_p3.debug:" >> $name.log
    cat ${name}_p3.debug >> $name.log
    echo "Contents of ${name}_p4.debug:" >> $name.log
    cat ${name}_p4.debug >> $name.log
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements, test skipped"
    exit 77
fi
