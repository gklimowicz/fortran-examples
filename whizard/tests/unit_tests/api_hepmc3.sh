#!/bin/sh
### Check WHIZARD command setup
echo "Running script $0"
if test -f HEPMC3_FLAG; then
    exec ./run_whizard_ut.sh --check api_hepmc
else
    echo "|=============================================================================|"
    echo "No HepMC3 available, test skipped"
    exit 77
fi
