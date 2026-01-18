#!/bin/sh
### Check WHIZARD module sf_base
echo "Running script $0"
if test -f LHAPDF5_FLAG; then
    exec ./run_whizard_ut.sh --check sf_lhapdf
else
    echo "|=============================================================================|"
    echo "No LHAPDF v4 or v5 available, test skipped"
    exit 77
fi
