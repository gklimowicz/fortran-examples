#!/bin/sh
echo "Running script $0"
name=`basename @script@`
if test -f ref-output/$name.ref; then
  if test -f OCAML_FLAG -a -f OPENLOOPS_FLAG -a -f FASTJET_FLAG -a -f LHAPDF6_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    mv $name.log $name.log.tmp
    cat $name.log.tmp | sed -e 's/Loading library:.*/Loading library: [...]/' > $name.log
    rm $name.log.tmp
    diff ref-output/$name.ref $name.log
    diffrc=$?
    grep --quiet "The desired process has not been found" $name.run.log
    greprc=$?
    grep --quiet "Requested library not installed" $name.run.log
    grep2rc=$?
    if test $diffrc -gt 0 -a $greprc -eq 0; then
      echo "|=============================================================================|"
      echo "OpenLoops process library missing"
      exit 99
    fi
    if test $diffrc -gt 0 -a $grep2rc -eq 0; then
      echo "|=============================================================================|"
      echo "OpenLoops process library missing"
      exit 99
    fi
    grep --quiet "'LUXqed_plus_PDF4LHC15_nnlo_100' not found" $name.run.log
    greprc=$?
    if test $diffrc -gt 0 -a $greprc -eq 0; then
      echo "|=============================================================================|"
      echo "LHAPDF: Data file 'LUXqed_plus_PDF4LHC15_nnlo_100' missing"
      exit 99
    fi
    exit $diffrc
  else
    echo "|=============================================================================|"
    echo "No O'Mega/OpenLoops matrix elements / FastJet / LHAPDF6 available"
    exit 77
  fi
else
  echo "|=============================================================================|"
  echo "$name.ref not found"
  exit 77
fi

