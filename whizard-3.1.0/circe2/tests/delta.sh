#! /bin/sh
########################################################################
if test -f OCAML_FLAG; then
    exec ./test_wrapper.sh @name@ check_generation </dev/null >/dev/null 2>&1
else
    echo "|=============================================================================|"
    echo "No OCaml for testing tools available, test skipped"
    exit 77
fi
