#!/bin/sh
### Check an LO RECOLA process with at least one color flow
echo "Running script $0"
if test -f OCAML_FLAG -a -f RECOLA_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    mv $name.log $name.log.tmp
    cat $name.log.tmp | sed -e 's/Loading library:.*/Loading library: [...]/' > $name.log
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega or RECOLA matrix elements available, test skipped"
    exit 77
fi
