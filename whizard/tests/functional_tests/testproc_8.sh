#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
./run_whizard.sh @script@ --no-logging
script=`basename @script@`
echo "Contents of ${script}_p.lha:" >> $script.log
cat ${script}_p.lha >> $script.log
diff ref-output/$script.ref $script.log
