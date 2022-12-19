#!/bin/sh
### Check WHIZARD command setup
echo "Running script $0"
if test -f LCIO_FLAG; then
    exec ./run_whizard_ut.sh --check api_lcio
else
    echo "|=============================================================================|"
    echo "No LCIO available, test skipped"
    exit 77
fi
