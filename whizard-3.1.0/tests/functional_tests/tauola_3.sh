#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f PYTHIA6_FLAG; then
    rm -f @script@_lib.* @script@_p?.*
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    echo "STDHEP EV4 version:" >> $script.log
    ./stdhep_rd ${script}_p1.ev4.hep 1 >> ${script}.log
    mv ${script}.log ${script}.log.tmp
    cat ${script}.log.tmp | sed -e 's/total blocks.*/total blocks: [...]/' -e 's/WHIZARD 3.*/WHIZARD [version]/' -e 's/date: .*/date: [...]/' -e 's/[ -][0-9].[0-9][0-9]e-0[5-9]/ 0.00e+00/g' -e 's/[ -][0-9].[0-9][0-9]e-1[0-9]/ 0.00e+00/g' > ${script}.log
    diff ref-output/${script}.ref ${script}.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available and/or PYTHIA6 disabled, test skipped"
    exit 77
fi
