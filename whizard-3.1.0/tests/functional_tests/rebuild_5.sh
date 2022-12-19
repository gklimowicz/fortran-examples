#!/bin/sh

echo "Running script $0"
script=`basename @script@`
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging --logfile=$script.1.log
    mv $script.dat $script.1.dat
    ./run_whizard.sh @script@ --no-logging --logfile=$script.2.log --no-rebuild
    mv $script.dat $script.2.dat
    diff $script.1.dat $script.2.dat
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
