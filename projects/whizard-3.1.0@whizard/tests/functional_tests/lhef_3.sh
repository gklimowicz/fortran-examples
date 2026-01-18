#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    rm -f @script@_lib.* @script@_p?.*
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model QCD
    echo "LHEF file contents (no beams):" >> $s.log
    cat ${s}.1.lhe | sed -e 's/^  <generator_version>.*$/[...]/' >> $s.log
    echo "LHEF file contents (with beams):" >> $s.log
    cat ${s}.2.lhe | sed -e 's/^  <generator_version>.*$/[...]/' >> $s.log
    echo "LHEF file contents (with beams and remnants as real):" >> $s.log
    cat ${s}.3.lhe | sed -e 's/^  <generator_version>.*$/[...]/' >> $s.log
    diff ref-output/$s.ref $s.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    

