#!/bin/sh

echo "Running script $0"
script=`basename @script@`
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging --logfile=$script.1.log
    ./run_whizard.sh @script@ --no-logging --logfile=$script.2.log --no-rebuild
    ./run_whizard.sh @script@ --no-logging --logfile=$script.3.log --no-rebuild --recompile
    cp ${script}p_i1.lo ${script}p_i1.tmp.lo
    touch ${script}p_i1.f90
    ./run_whizard.sh @script@ --no-logging --logfile=$script.4.log --no-rebuild --recompile
    cat $script.1.log $script.2.log $script.3.log $script.4.log > $script.log
    newerfile=$(ls -t ${script}p_i1.lo ${script}p_i1.tmp.lo | head -n1)
    if test "$newerfile" = ${script}p_i1.lo; then
	echo "Library was updated" >> $script.log
    else
	echo "Library was not updated" >> $script.log
    fi
    rm -f ${script}p_i1.tmp.lo
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
