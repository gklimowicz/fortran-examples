#!/bin/sh
### Check WHIZARD eio_hepmc module setup
echo "Running script $0"
if test -f HEPMC3_FLAG; then
    exec ./run_whizard_ut.sh --check eio_hepmc
else
    echo "|=============================================================================|"
    echo "No HepMC3 available, test skipped"
    exit 77
fi
