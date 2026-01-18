#!/bin/sh
### Check WHIZARDs final state shower/matching
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    script=`basename @script@`
    echo "Contents of ${script}_1.evt:" >> $script.log
    cat ${script}_1.evt >> $script.log
    echo "Contents of ${script}_2.evt:" >> $script.log
    cat ${script}_2.evt >> $script.log
    diff ref-output/$script.ref $script.log

else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi

