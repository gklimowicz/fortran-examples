#!/bin/sh
### Check WHIZARD hepmc module setup
echo "Running script $0"
if test -f HEPMC3_FLAG; then
    exec ./run_whizard_ut.sh --check hepmc
else
    echo "|=============================================================================|"
    echo "No HepMC3 available, test skipped"
    exit 77
fi
