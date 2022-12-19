#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model QCD
    echo "Contents of ${name}a.mokka.evt" >> $name.log
    cat ${name}a.mokka.evt >> $name.log
    echo "Contents of ${name}b.mokka.evt" >> $name.log
    cat ${name}b.mokka.evt >> $name.log
    cat $name.log | sed -e 's/WHIZARD 3.*$/WHIZARD [version]/' > $name.log.tmp
    mv $name.log.tmp $name.log
    diff ref-output/`basename @script@`.ref `basename @script@`.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
