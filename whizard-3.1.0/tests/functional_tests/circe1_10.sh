#!/bin/sh
### Check WHIZARD/O'Mega Circe/ISR setup
echo "Running script $0"
if test -f OCAML_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model SM
    echo "Contents of ${script}_p1_1.evt" >> $script.log
    cat ${script}_p1_1.evt >> $script.log
    echo "Contents of ${script}_p1_2.evt" >> $script.log
    cat ${script}_p1_2.evt >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
