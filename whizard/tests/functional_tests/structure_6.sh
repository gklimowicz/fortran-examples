#!/bin/sh
### Check WHIZARD control structures
echo "Running script $0"
name=`basename @script@`
./run_whizard.sh @script@ --no-logging --no-model
echo "Contents of ${name}a.out:" >> $name.log
cat ${name}a.out >> $name.log
echo "Contents of ${name}b.out:" >> $name.log
cat ${name}b.out >> $name.log
diff ref-output/$name.ref $name.log
