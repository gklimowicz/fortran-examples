#!/bin/sh
### Check WHIZARD beam setup
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging --no-model
    script=`basename @script@`
    for i in 1 2 3 4; do
	echo "Contents of ${script}p_$i.mokka.evt:" >> $script.log
	cat ${script}p_$i.mokka.evt >> $script.log
    done
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
