#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f LCIO_FLAG; then
    rm -f @script@_lib.* @script@_p?.*
    s=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model QED
    echo "Output from running ${s}_rd:" >> ${s}.log
    ./lcio_rd ${s}_p.slcio 21 1 >> ${s}.log
    cat ${s}.log | sed -e 's/^ date:.*$/ date: [...]/' | sed -e 's/timestamp .*$/ timestamp [...]/' > ${s}.log.tmp && mv ${s}.log.tmp ${s}.log
    diff ref-output/$s.ref ${s}.log    
else
    echo "|=============================================================================|"
    echo "No LCIO or no O'Mega matrix elements available, test skipped"
    exit 77
fi

