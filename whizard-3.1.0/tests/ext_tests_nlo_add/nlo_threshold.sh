#!/bin/sh
# time ~
echo "Running script $0"
if test -f OCAML_FLAG -a -f OPENLOOPS_FLAG; then
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --debug threshold
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
  else
    echo "|=============================================================================|"
    echo "No O'Mega and/or OpenLoops matrix elements available"
    exit 77
fi
