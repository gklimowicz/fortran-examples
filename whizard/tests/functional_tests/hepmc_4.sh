#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f HEPMC2_FLAG || test -f OCAML_FLAG -a -f HEPMC3_FLAG; then
    rm -f @script@_lib.* @script@_p?.*
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model QED
    # HepMC file contents depend on the specific of the installation
    # for i in 0 1 2; do
    # 	echo >> $s.log
    # 	echo "HepMC file contents:" >> $s.log
    # 	cat ${s}_p.$i.hepmc >> $s.log
    # done
    diff ref-output/$s.ref $s.log
else
    echo "|=============================================================================|"
    echo "No HepMC available, test skipped"
    exit 77
fi

