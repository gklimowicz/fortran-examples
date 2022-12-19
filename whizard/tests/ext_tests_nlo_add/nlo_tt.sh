#!/bin/bash
# time ~9min
echo "Running script $0"
set -o pipefail # set return code of a pipe to last non-zero exit code
if test -f OCAML_FLAG -a -f OPENLOOPS_FLAG -a -f HEPMC2_FLAG; then
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
    @python_bin@ @share_dir@/compare-integrals.py $s \
      @share_dir@/extra_integration_results.dat | tee -a $s.run.log
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
    @python_bin@ @share_dir@/check-hepmc-weights.py $s | tee -a $s.run.log
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
  else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available"
    exit 77
fi
