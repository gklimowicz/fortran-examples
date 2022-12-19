#!/bin/sh
# time ~
echo "Running script $0"
if test -f OCAML_FLAG -a -f OPENLOOPS_FLAG; then
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
    mv $s.log $s.log.tmp
    cat $s.log.tmp | sed -e 's/Loading library:.*/Loading library: [...]/' > $s.log
    diff ref-output/$s.ref $s.log
  else
    echo "|=============================================================================|"
    echo "No O'Mega and/or OpenLoops matrix elements available"
    exit 77
fi
