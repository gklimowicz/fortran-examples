#!/bin/sh
### Testing the integration of the virtual component of pp -> Zj
### and the simulation of events.
echo "Running script $0"
if test -f OCAML_FLAG -a -f OPENLOOPS_FLAG -a -f FASTJET_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    mv $name.log $name.log.tmp
    cat $name.log.tmp | sed -e 's/Loading library:.*/Loading library: [...]/' > $name.log
    cat ${name}_p1_fks_regions.out >> $name.log
    echo "Contents of ${name}_p1.debug:" >> $name.log
    cat ${name}_p1.debug >> $name.log
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega/OpenLoops matrix elements / FastJet available. Test skipped."
    exit 77
fi
