#!/bin/sh
### Check WHIZARD for a simple test process
script=`basename @script@`
echo "Running script $0"
./run_whizard.sh @script@ --no-model
out0=`grep "Events: writing to raw file" ${script}.log`
if test -z "$out0"; then
    exit 1
fi
./run_whizard.sh @script@ --no-model --no-rebuild
out1=`grep "Events: reading from raw file" ${script}.log`
out2=`grep "Events: parameter mismatch" ${script}.log`
if test -n "$out1" -a -z "$out2"; then
    exit 0
else
    exit 2
fi
