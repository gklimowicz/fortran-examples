#!/bin/sh
### Check WHIZARD/O'Mega interaction with SUSYHIT
echo "Running script $0"
PRG=susyhit
if (which $PRG >/dev/null 2>&1); then
    name=`basename @script@`
    ./run_whizard.sh @script@ --no-logging
    cat ${name}.log | sed -e 's/| SLHA: [23].*/| SLHA: [version]/' > $name.log.tmp
    mv $name.log.tmp $name.log
    if test -f suspect2.out -a -f slhaspectrum.in -a -f susyhit_slha.out; then
	diff ref-output/`basename @script@`.ref `basename @script@`.log
    else
	exit 1
    fi
else
  echo "|=============================================================================|"
  echo "$PRG executable not found in PATH, test skipped"
  exit 77
fi
