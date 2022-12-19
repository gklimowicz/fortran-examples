#!/bin/sh
### Check WHIZARD/O'Mega: ISR/EPA handler
echo "Running script $0"
if test -f OCAML_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model SM
    ./`basename @script@`_count
    cat `basename @script@`_count.out >> `basename @script@`.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
