#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f HEPMC2_FLAG || test -f OCAML_FLAG -a -f HEPMC3_FLAG; then
    rm -f @script@_lib.* @script@_p?.*
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model QCD
    cat $name.log | sed -e 's/WHIZARD 2.*$/WHIZARD [version]/' > $name.log.tmp
    mv $name.log.tmp $name.log
    echo "Contents of ${name}a.weights.dat" >> $name.log
    cat ${name}a.weights.dat >> $name.log
    echo "Contents of ${name}b.weights.dat" >> $name.log
    cat ${name}b.weights.dat >> $name.log
    diff ref-output/`basename @script@`.ref `basename @script@`.log
else
    echo "|=============================================================================|"
    echo "HepMC disabled or no O'Mega matrix elements available, test skipped"
    exit 77
fi
    
