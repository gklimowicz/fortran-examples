#!/bin/sh
### Check WHIZARD POWHEG matching with dummy virtual matrix-elements
echo "Running script $0"
if test -f OCAML_FLAG; then
    name=`basename @script@`
    rm -f ${name}_p1.pg
    rm -f ${name}_p2.pg
    ./run_whizard.sh @script@ --no-logging
    echo "Contents of ${name}_p1.debug:" >> $name.log
    cat ${name}_p1.debug >> $name.log
    echo "Contents of ${name}_p2.debug:" >> $name.log
    cat ${name}_p2.debug >> $name.log
    echo "Contents of ${name}_p1.pg:" >> $name.log
    cat ${name}_p1.pg >> $name.log
    echo "Contents of ${name}_p2.pg:" >> $name.log
    cat ${name}_p2.pg >> $name.log
    diff -b ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
