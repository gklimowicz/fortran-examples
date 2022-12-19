#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f LCIO_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    echo "Contents of ${name}a.evt" >> $name.log
    cat ${name}a.evt | grep 'prt(' \
	| sed -e 's/^\(.*|.*|\).*\(|.*\)$/\1 *** \2/' >> $name.log
    echo "Contents of ${name}b.evt" >> $name.log
    cat ${name}b.evt | grep 'prt(' \
	| sed -e 's/^\(.*|.*|\).*\(|.*\)$/\1 *** \2/' >> $name.log
    diff ref-output/`basename @script@`.ref `basename @script@`.log
else
    echo "|=============================================================================|"
    echo "No LCIO or no O'Mega matrix elements available, test skipped"
    exit 77
fi
    
