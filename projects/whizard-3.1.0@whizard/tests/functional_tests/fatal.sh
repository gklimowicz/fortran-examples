#!/bin/sh
### Check WHIZARD error and fatal error handling. Error and
### fatal error output is compared 
echo "Running script $0"
./run_whizard.sh @script@ --no-logging 
if test $? -eq 1; then
    diff ref-output/`basename @script@`.ref `basename @script@`.log
else
    exit 42
fi


