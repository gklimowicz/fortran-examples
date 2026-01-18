#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    rm -f ${script}_lib*
    rm -f ${script}_x*
    ./run_whizard.sh @script@ --no-logging --no-model -J "4711"
    script=`basename @script@`
    echo "* Files created by compile:" >> $script.log
    ls ${script}_lib.f90 >> $script.log
    ls ${script}_lib.la >> $script.log
    ls ${script}_lib.makefile >> $script.log
    ls ${script}_x_i1.f90 >> $script.log
    ls ${script}_x_i1.lo >> $script.log
    echo "* Files created by integrate:" >> $script.log
    ls ${script}_x.*.1.* >> $script.log
    echo "* Files created by simulate:" >> $script.log
    ls ${script}_x.*.2.* >> $script.log
    echo "* Files created by analysis:" >> $script.log
    ls ${script}_x.*.3.* >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
