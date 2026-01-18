#!/bin/sh
### Check WHIZARD unpack/pack feature
echo "Running script $0"

if test -f GZIP_FLAG; then
    script=`basename @script@`
    tar_create="tar -czf"
    tar_list="tar -tzf"

    rm -rf $script.1 $script.2 $script.3 $script.4
    rm -rf $script.1.tgz $script.2.tgz $script.3.tgz $script.4.tgz
    
    mkdir $script.1
    touch $script.1/foo
    $tar_create $script.1.tgz $script.1
    rm -rf $script.1
    
    mkdir $script.2
    touch $script.2/bar
    $tar_create $script.2.tgz $script.2
    rm -rf $script.2
    
    mkdir $script.3
    touch $script.3/foo

    mkdir $script.4
    touch $script.4/bar
    
    ./run_whizard.sh @script@ --no-model --no-logging --no-library \
	--unpack $script.1.tgz --pack $script.3 \
	--unpack $script.2.tgz --pack $script.4
    
    echo "Contents of directory $script.1:" > $script.log
    ls $script.1 >> $script.log

    echo "Contents of directory $script.2:" >> $script.log
    ls $script.2 >> $script.log
    
    echo "Contents of file $script.3.tgz:" >> $script.log
    $tar_list $script.3.tgz >> $script.log
    
    echo "Contents of file $script.4.tgz:" >> $script.log
    $tar_list $script.4.tgz >> $script.log
    
    diff ref-output/$script.ref $script.log
else
    echo "|=============================================================================|"
    echo "gzip unavailable, test skipped"
    exit 77
fi

