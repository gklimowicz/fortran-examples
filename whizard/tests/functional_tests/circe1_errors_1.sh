#!/bin/sh
### Check Circe1 Error handling
echo "Running script $0"
if test -f OCAML_FLAG; then
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model SM \
	| grep '^ *circe1:' | grep -v '\$Id:' >${script}_circe1.log
    diff ref-output/$script.ref ${script}_circe1.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
