#!/bin/sh
### Check WHIZARD for a simple test process
name=`basename @script@`
echo "Running script $0"
if test -f EVENT_ANALYSIS_FLAG; then
    ./run_whizard.sh @script@ --no-logging --no-model
    echo "Plot range settings in ${name}_plots.tex:" >> $name.log
    cat ${name}_plots.tex | grep graphrange >> $name.log
    diff ref-output/`basename @script@`.ref `basename @script@`.log
else
    echo "|=============================================================================|"
    echo "No event analysis PS and PDF files can be generated, test skipped"
    exit 77
fi
    
