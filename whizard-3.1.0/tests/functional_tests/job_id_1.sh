#!/bin/sh
### Check WHIZARD with Job ID
echo "Running script $0"
./run_whizard.sh @script@ --no-logging --no-model --job-id="JOB_ID_STR"
diff ref-output/`basename @script@`.ref `basename @script@`.log
