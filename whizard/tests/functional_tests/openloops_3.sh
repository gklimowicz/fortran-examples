#!/bin/sh
### Check WHIZARD NLO Event Generation from separate components
echo "Running script $0"
if test -f OCAML_FLAG -a -f OPENLOOPS_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    cat ${name}_p1.debug >> $name.log
    cat ${name}_p2.debug >> $name.log
    cat ${name}_p3.debug >> $name.log
    mv $name.log $name.log.tmp
    cat $name.log.tmp | sed -e 's/Loading library:.*/Loading library: [...]/' > $name.log
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
