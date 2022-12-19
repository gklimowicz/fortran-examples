#!/bin/sh
# time ~
echo "Running script $0"
if test -f OCAML_FLAG -a -f OPENLOOPS_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --debug threshold
    whizardrc=$?
    grep --quiet "\[OpenLoops\] Requested library not installed." $name.run.log
    greprc=$?
    if test $whizardrc -gt 0 -a $greprc -eq 0; then
      echo "|=============================================================================|"
      echo "OpenLoops process library missing"
      exit 77
    fi
    exit $whizardrc
  else
    echo "|=============================================================================|"
    echo "No O'Mega and/or OpenLoops matrix elements available"
    exit 77
fi
