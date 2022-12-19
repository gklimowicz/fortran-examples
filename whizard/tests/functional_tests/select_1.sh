#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    script=`basename @script@`
    mv $script.log $script.log.tmp
    cat $script.log.tmp | sed -e 's/WHIZARD 3.*$/WHIZARD [version]/' > $script.log
    rm -f $script.log.tmp
    echo "Contents of ${script}a.dat:" >> $script.log
    cat ${script}a.dat >> $script.log
    echo "Contents of ${script}b.dat:" >> $script.log
    cat ${script}b.dat >> $script.log
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
