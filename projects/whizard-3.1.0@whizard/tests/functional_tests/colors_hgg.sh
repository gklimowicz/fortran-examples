#!/bin/sh
### Check WHIZARD/O'Mega with color correlations
echo "Running script $0"
if test -f OCAML_FLAG; then
    ./run_whizard.sh @script@ --no-logging
    if test $? -eq 5; then
	diff ref-output/`basename @script@`.ref `basename @script@`.log
    fi
else
    echo "|=============================================================================|"
    echo "No O'Mega matrix elements available, test skipped"
    exit 77
fi
