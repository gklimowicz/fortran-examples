#!/bin/sh
### Check WHIZARD command setup
echo "Running script $0"
if test -f HEPMC2_FLAG; then
    exec ./run_whizard_ut.sh --check api_hepmc
else
    echo "|=============================================================================|"
    echo "No HepMC2 available, test skipped"
    exit 77
fi
