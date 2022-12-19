#!/bin/sh
### Check WHIZARD compilations module
echo "Running script $0"
if test -f STATIC_FLAG; then
    if test -f OCAML_FLAG; then
	exec ./run_whizard_ut.sh --check compilations_static
    else
	echo "|=============================================================================|"
	echo "No O'Mega matrix elements available, test skipped"
	exit 77
    fi
else
    echo "|=============================================================================|"
    echo "Static libraries disabled, test skipped"
    exit 77
fi
