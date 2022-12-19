#!/bin/bash
### download script for external data needed for 'SM_tt_threshold' model
### when coninuum matching is enabled

DATFILE="SM_tt_threshold_Vmatrices.tar.gz"

wget http://whizard.hepforge.org/versions/data_files/$DATFILE
tar -xvf $DATFILE

if [ "`pwd`" != "`readlink -f $BASH_SOURCE | xargs dirname`" ]; then
  mv ${DATFILE%.tar.gz} `dirname $BASH_SOURCE`
fi

rm -f $DATFILE
