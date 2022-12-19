#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    rm -f @script@_lib.* @script@_p?.*
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model QCD
    for i in 0 1 2; do
	echo >> $s.log
	echo "LHEF file contents:" >> $s.log
	cat ${s}_p.$i.lhe | sed -e 's/^<generator version=.*$/<generator version=[...]>WHIZARD<\/generator>/' >> $s.log
    done
    diff ref-output/$s.ref $s.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
