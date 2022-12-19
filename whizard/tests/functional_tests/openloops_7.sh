#!/bin/sh
### Check a beam-polarized cross section computation with OpenLoops
echo "Running script $0"
if test -f OCAML_FLAG -a -f OPENLOOPS_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    mv $name.log $name.log.tmp
    cat $name.log.tmp | sed -e 's/Loading library:.*/Loading library: [...]/' > $name.log
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega or OpenLoops matrix elements available, test skipped"
    exit 77
fi
