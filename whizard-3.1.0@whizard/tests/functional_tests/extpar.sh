#!/bin/sh
### Check WHIZARD/O'Mega with some simple QED process
echo "Running script $0"
./run_whizard.sh @script@ --no-logging --model Test
diff ref-output/`basename @script@`.ref `basename @script@`.log
