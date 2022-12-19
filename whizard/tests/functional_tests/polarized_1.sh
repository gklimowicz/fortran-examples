#!/bin/sh
echo "Running script $0"
if test -f OCAML_FLAG; then
    rm -f @script@_lib.* @script@_p?.*
    script=`basename @script@`    
    ./run_whizard.sh @script@ --no-logging --model SM
    echo "Contents of ${script}_p1.hepevt:" >> ${script}.log
    cat ${script}_p1.hepevt >> ${script}.log
    echo "STDHEP EV4 version:" >> $script.log
    ./stdhep_rd ${script}_p1.ev4.hep 1 >> ${script}.log
    mv ${script}.log ${script}.log.tmp    
    cat ${script}.log.tmp | sed -e 's/total blocks.*/total blocks: [...]/' -e 's/WHIZARD 3.*/WHIZARD [version]/' -e 's/date: .*/date: [...]/' > ${script}.log
    diff ref-output/${script}.ref ${script}.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
