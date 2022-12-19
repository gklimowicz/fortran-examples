#!/bin/sh
### Check WHIZARD decays module
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard_ut.sh --check shower
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi

