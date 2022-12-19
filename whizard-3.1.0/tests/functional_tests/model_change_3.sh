#!/bin/sh
### Check WHIZARD model setup
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    name=`basename @script@`
    grep '| There were' $name.log > $name.tmp
    mv $name.tmp $name.log
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi

