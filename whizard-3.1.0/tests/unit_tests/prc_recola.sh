#!/bin/sh
### Check WHIZARD module RECOLA
echo "Running script $0"
if test -f OCAML_FLAG -a -f RECOLA_FLAG; then
    exec ./run_whizard_ut.sh --check prc_recola
else
    echo "|=============================================================================|"
    echo "No RECOLA available, test skipped"
    exit 77
fi
