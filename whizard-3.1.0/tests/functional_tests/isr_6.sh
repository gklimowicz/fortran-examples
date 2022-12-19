#!/bin/sh
### Check WHIZARD/O'Mega with externally linked LHAPDF library
echo "Running script $0"
if test -f OCAML_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model SM
    ./${script}_digest
    echo "Contents of file ${script}_digest.out:" >> $script.log
    cat ${script}_digest.out >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
