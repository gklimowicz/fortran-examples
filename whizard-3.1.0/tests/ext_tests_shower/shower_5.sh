#!/bin/bash
# time ? min
echo "Running script $0"
set -o pipefail # set return code of a pipe to last non-zero exit code
if test -f OCAML_FLAG -a -f PYTHIA6_FLAG -a -f PYTHON_FLAG; then
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
    @python_bin@ @share_dir@/check-debug-output.py ${s}_uu | tee -a $s.run.log
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
    @python_bin@ @share_dir@/check-debug-output.py ${s}_dd | tee -a $s.run.log
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
    @python_bin@ @share_dir@/check-debug-output.py ${s}_cc | tee -a $s.run.log
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
    @python_bin@ @share_dir@/check-debug-output.py ${s}_ss | tee -a $s.run.log
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
    @python_bin@ @share_dir@/check-debug-output.py ${s}_bb | tee -a $s.run.log
    rc=$?; if [ $rc != 0 ]; then exit $rc; fi
  else
    echo "|=============================================================================|"
    echo "O'Mega, PYTHIA6 and/or Python disabled/not available, test skipped"
    exit 77
fi