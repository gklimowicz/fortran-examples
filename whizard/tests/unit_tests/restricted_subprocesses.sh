#!/bin/sh
### Check WHIZARD restricted_subprocesses module
echo "Running script $0"
if test -f OCAML_FLAG; then
    exec ./run_whizard_ut.sh --check restricted_subprocesses
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
