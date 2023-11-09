#!/bin/sh
### Check WHIZARD self-built static executable
echo "Running script $0"
script=`basename @script@`
exe=$script.exe
if test -f STATIC_FLAG; then
    if test -f OCAML_FLAG; then
	rm -f $exe
	./run_whizard.sh @script@ --no-logging --no-model 
	if test -x "$exe"; then
	    ./$exe $exe.sin --logfile $exe.log --no-library --no-logging --no-model --no-banner --rebuild
	    echo >> $script.log
	    cat $exe.log >> $script.log
	fi
	diff ref-output/$script.ref $script.log
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
