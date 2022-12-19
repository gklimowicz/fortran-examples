#!/bin/sh
### Check WHIZARD control structures
echo "Running script $0"
script=`basename @script@`
./run_whizard.sh @script@ \
    -e "'int i = 1'" \
    -f ${script}_a.sin \
    --execute \'i = 4  a = 12\' \
    --file ${script}_b.sin \
    -e \"int q = 3\" \
    --no-logging --no-model --no-library
diff ref-output/$script.ref $script.log
