#!/bin/sh
### Check WHIZARD process library setup
echo "Running script $0"
if test -f OCAML_FLAG -a -f EVENT_ANALYSIS_FLAG; then
    exec ./run_whizard_ut.sh --check prc_omega_diags
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, and or no event analysis"
    echo "     PS and PDF files can be generated, test skipped"
    exit 77
fi
    
