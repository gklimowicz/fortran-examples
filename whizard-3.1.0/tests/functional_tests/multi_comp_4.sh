#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging --no-model
    script=`basename @script@`
    evtfile="${script}p.short.evt"
    echo "Number of mu+ events  = " \
	`grep "^ -13" $evtfile | wc | awk '{print $1}'` >> $script.log 
    echo "Number of tau+ events = " \
	`grep "^ -15" $evtfile | wc | awk '{print $1}'` >> $script.log 
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi    
