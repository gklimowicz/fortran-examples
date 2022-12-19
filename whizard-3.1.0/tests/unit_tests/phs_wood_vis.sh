#!/bin/sh
### Check WHIZARD integrations history files
echo "Running script $0"
if test -f EVENT_ANALYSIS_FLAG; then
    exec ./run_whizard_ut.sh --check phs_wood_vis
else
    echo "|=============================================================================|"
    echo "No event analysis PS and PDF files can be generated, test skipped"
    exit 77
fi
    

