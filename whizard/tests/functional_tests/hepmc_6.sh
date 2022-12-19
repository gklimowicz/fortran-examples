#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f HEPMC2_FLAG || test -f OCAML_FLAG -a -f HEPMC3_FLAG; then
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model QED
    echo "Output from running ${s}_rd:" > ${s}.log
    ./${s}_rd >> ${s}.log
    diff ref-output/$s.ref ${s}.log
else
    echo "|=============================================================================|"
    echo "No HepMC available, test skipped"
    exit 77
fi

