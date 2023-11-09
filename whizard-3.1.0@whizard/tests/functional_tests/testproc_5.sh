#!/bin/sh
### Check WHIZARD for a simple test process
name=`basename @script@`
echo "Running script $0"
./run_whizard.sh @script@ --no-logging
echo "Contents of ${name}a.weights.dat" >> $name.log
cat ${name}a.weights.dat >> $name.log
echo "Contents of ${name}b.weights.dat" >> $name.log
cat ${name}b.weights.dat >> $name.log
echo "Contents of ${name}c.weights.dat" >> $name.log
cat ${name}c.weights.dat >> $name.log
diff ref-output/$name.ref $name.log
