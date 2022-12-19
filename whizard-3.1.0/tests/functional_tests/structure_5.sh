#!/bin/sh
### Check WHIZARD control structures
echo "Running script $0"
./run_whizard.sh @script@ --no-logging --no-model
diff ref-output/`basename @script@`.ref `basename @script@`.log
