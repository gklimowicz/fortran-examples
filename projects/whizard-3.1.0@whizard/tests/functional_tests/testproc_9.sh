#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
./run_whizard.sh @script@ --no-logging
script=`basename @script@`
for i in 1 2 3; do
    echo >> $script.log
    echo "Contents of ${script}_$i.debug:" >> $script.log
    cat ${script}_$i.debug >> $script.log
done
diff ref-output/$script.ref $script.log
