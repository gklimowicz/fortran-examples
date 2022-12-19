#!/bin/sh
### Check WHIZARD/O'Mega Circe/EPA setup
echo "Running script $0"
if test -f OCAML_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model SM
    ### We don't know the valid output yet
    # if test $? -eq 1; then
    # 	diff ref-output/$script.ref $script.log
    # else
    # 	exit 42
    # fi
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
