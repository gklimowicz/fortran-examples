#!/bin/sh
### Check WHIZARD for a simple test process
echo "Running script $0"
if test -f OCAML_FLAG -a -f HEPMC2_FLAG || test -f OCAML_FLAG -a -f HEPMC3_FLAG; then
    rm -f @script@_lib.* @script@_p?.*
    script=`basename @script@`
    ./run_whizard.sh @script@ --no-logging --model QCD
    # The HEPEVT output, extract the W helicity for each event
    echo "W polarization in ${script}_p.hepevt:" >> $script.log
    grep ' -\?24 ' ${script}_p.hepevt | sed -e 's/^P .\+ \(-\?24\) .\+ \(-\?[01]\)$/\1 \2/' >> $script.log
    # The HepMC output may be version-dependent, so don't check automatically.
    #echo "W polarization in ${script}_p.hepmc:" >> $script.log
    #grep ' -\?24 ' ${script}_p.hepmc >> $script.log
    diff ref-output/${script}.ref $script.log
else
    echo "|=============================================================================|"
    echo "HepMC disabled or no O'Mega matrix elements available, test skipped"
    exit 77
fi
    
