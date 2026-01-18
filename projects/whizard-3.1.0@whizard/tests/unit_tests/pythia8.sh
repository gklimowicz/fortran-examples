#!/bin/sh
### Check WHIZARD module pythia8
echo "Running script $0"
if test -f PYTHIA8_FLAG; then
    exec ./run_whizard_ut.sh --check pythia8
else
    echo "|=============================================================================|"
    echo "No Pythia8 available, test skipped"
    exit 77
fi
