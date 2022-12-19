#!/bin/sh
# Testing event generation with a multi-flavor final state
# at next-to-leading order QCD
echo "Running script $0"
if test -f OCAML_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    echo "Contents of ${name}_p1.debug:" >> $name.log
    cat ${name}_p1.debug >> $name.log
    cat ${name}_p1_fks_regions.out >> $name.log
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
