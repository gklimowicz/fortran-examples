#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging --no-model
    script=`basename @script@`
    evtfile="${script}_p1.short.evt"
    echo "Number of mu+,mu- events     = " \
	`grep "^ 13" $evtfile | wc | awk '{print $1}'` >> $script.log 
    echo "Number of tau+,tau- events   = " \
	`grep "^ 15" $evtfile | wc | awk '{print $1}'` >> $script.log 
    echo "Number of gamma,gamma events = " \
	`grep "^ 22" $evtfile | wc | awk '{print ($1 / 2)}'` >> $script.log 
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
    
