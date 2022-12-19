#!/bin/sh
### Testing the integration of the real component of ee -> jjj
### and the simulation of events.
echo "Running script $0"
if test -f OCAML_FLAG -a -f OPENLOOPS_FLAG -a -f FASTJET_FLAG; then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    mv $name.log $name.log.tmp
    cat $name.log.tmp | sed -e 's/Loading library:.*/Loading library: [...]/' > $name.log
    cat ${name}_p1_fks_regions.out >> $name.log
    echo "Contents of ${name}_p1.debug:" >> $name.log
    cat ${name}_p1.debug | sed -e 's/\(sqme_rad =  [0-9]\.[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\)[0-9]E-\([0-9][0-9]\)/\1XE-\2/' -e '/prt(o:21/s/[0-9]E/XE/g' >> $name.log
    diff ref-output/$name.ref $name.log
else
    echo "|=============================================================================|"
    echo "No O'Mega/OpenLoops matrix elements / FastJet available. Test skipped."
    exit 77
fi
