#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
./run_whizard.sh @script@ --no-logging
diff ref-output/`basename @script@`.ref `basename @script@`.log
    
