#!/usr/bin/perl
## setup_e - set up GCM modelE run
## This script is not supposed to be executed by users but rather
## is a part of Makefile functionality. It is run as a part of
##           gmake setup RUN=RunID
## Inside the Makefile it is called (from .../decks): setup_e RunID
## It looks for $HOME/.modelErc file and extracts all necessary options from
## there. If such file is not present the default options which are
## specified below will be used. Those default options are adjusted for
## working environment of ra.giss.nasa.gov.
##
## The options setup_e is looking for are:
## CMRUNDIR - directory to which all run directories will be linked.
## EXECDIR - path to directory with modelE scripts and with some executables.
## SAVEDISK - a directory (big) where all run directories will be created.
## GCMSEARCHPATH - directory to search for gcm input files.
## MAILTO - email address of the user. (default `whoami`)
## UMASK - the value of 'umask' to be used for model runs.

## default settings
$CMRUNDIR="/u/cmrun";
$SAVEDISK="/raid1";
$GCMSEARCHPATH="/u/cmrun";
$NETCDFHOME="/usr/bin";
$MAILTO="";
$UMASK=002;
$MIN_STACK=0;
$mpi=0;
$nproc=1;
$MPIDISTR="";
$MPIDIR="";
$MPIRUN_COMMAND="";
$LOCATION="";

## if $HOME/.modelErc is present get settings from there

if (exists $ENV{MODELERC}) {
    $modelerc = $ENV{MODELERC}
}
else {
    $modelerc = (getpwuid($>))[7]."/.modelErc";
}

if ( -f $modelerc ) {
    print "Using settings from $modelrc \n" ;
    open MODELERC, $modelerc or die "can't open $modelerc";
    while(<MODELERC>) {
	$CMRUNDIR = $1 if /^ *CMRUNDIR *= *(\S+)/;
	$SAVEDISK = $1 if /^ *SAVEDISK *= *(\S+)/;
	$GCMSEARCHPATH = $1 if /^ *GCMSEARCHPATH *= *(\S+)/;
	$MAILTO = $1 if /^ *MAILTO *= *(\S+)/;
	$UMASK = oct "$1" if /^ *UMASK *= *(\S+)/;
	$NETCDFHOME = $1 if /^ *NETCDFHOME *= *(\S+)/;
        $MPIDISTR = $1 if /^ *MPIDISTR *= *\"?([^ ^#][^#^"]*).*\n/;
        $MPIDIR = $1 if /^ *MPIDIR *= *\"?([^ ^#][^#^"]*).*\n/;
	$MPIRUN_COMMAND = $1 if /^ *MPIRUN_COMMAND *= *\"?([^ ^#][^#^"]*).*\n/;
	$LOCATION = $1 if /^ *LOCATION *= *\"?([^ ^#][^#^"]*).*\n/;
						#	])])]);
    }
    close MODELERC;
} else {
    print "$HOME/.modelErc is not present. Using default settings.\n";
}

if ( $MAILTO =~ /^\s*$/ ) { $MAILTO = `whoami`; chop $MAILTO; }
print "CMRUNDIR = $CMRUNDIR\n";
print "GCMSEARCHPATH = $GCMSEARCHPATH\n";
print "SAVEDISK = $SAVEDISK\n";
print "MAILTO = $MAILTO\n";
printf "UMASK = %03lo\n", $UMASK;
print "MPIRUN_COMMAND = $MPIRUN_COMMAND\n";

while ($_ = $ARGV[0], /^-/) {
    shift;
    last if /^--$/;
    if (/^-mpi\b/) { $mpi = 1; next;}
    if (/^-mpidistr\b/) { $MPIDISTR = shift; $mpi = 1; next;}
    print "setup: unknown option $_ \n"; exit 1;
}


if ( $#ARGV != 0 ) { 
    print "This script is not supposed to be run outside of Makefile\n";
    print "Inside the Makefile it is called:";
    print "setup_e.pl {-omp|-mpi} num_proc runID\n";
    print "$#ARGV\n";
    exit 1; 
}

$runID = shift;
$rfile = "$runID.R";
$RunDir = "$SAVEDISK/$runID";
$lockfile = "$RunDir/lock";
$CMRUN = $CMRUNDIR;
$DeckDir = `pwd`; chop $DeckDir;
$umask_inv = $UMASK ^ 0777;
$umask_str = sprintf "%lo", $UMASK;
$NETCDFBIN = "$NETCDFHOME/bin";
$netcdf_template_file = "$runID.nctemp";
if ( $MPIDIR ) { $MPIBIN = "$MPIDIR/bin/"; }
else { $MPIBIN = ""; }

## check if this run is already running
if ( -f $lockfile ) {
    print "            **********************                \n";
    print "$1 seems to be already running in $SAVEDISK/$runID\n";
    print "If you think it is an error, then most probably this\n";
    print "task was interrupted in an unusual way. Please check.\n";
    print "Then remove the lock file:\n";
    print "$lockfile\n";
    print "and restart the setup.\n";
    print "            **********************                \n";
    exit 1;
}

## check if the rundeck is present
if ( ! -s $rfile ) {
   print "File does not exist: $rfile\n";
   exit 1;
}

print "setting up run $runID\n";
print "output files will be saved in $SAVEDISK/$runID\n";

## Create the run directory if necessary
if ( ! -d $RunDir ) {
    mkdir $RunDir, 0777 & $umask_inv 
	or die "Can't create $RunDir. Aborting setup.\n";
}

## Check that link is not already correct (added by gavin)
if ( -l $runID ) {
    if ( `ls -l $runID` !~ /-> *$RunDir$/ ) {
	print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	print "A link in your local directory \n";
	print "    ./$runID\n";
	print "exists and is pointing to something else.\n";
	print "Will not overwrite it.\n";
	print "Please check. Aborting setup.\n";
	print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	exit 1;
    }
} elsif ( -e $runID ) {
	print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	print "An object in your local directory \n";
	print "    ./$runID\n";
	print "exists and is not a symbolic link.\n";
	print "Will not overwrite it.\n";
	print "Please check. Aborting setup.\n";
	print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	exit 1;
} else {
    symlink $RunDir, $runID or die "Can't create link $runID in local dir";
}

## Also link to $CMRUN if different from $SAVEDISK (added by gavin)
if ( $CMRUN ne $SAVEDISK ) {
    if ( -l "$CMRUN/$runID" ) {                 # is a link
	if (`ls -l $CMRUN/$runID` !~ /-> *$RunDir$/ ) {
	    unlink "$CMRUN/$runID" or die "Can't rm old link $CMRUN/$runID";
	    symlink "$RunDir", "$CMRUN/$runID" or die "Can't create link $CMRUN/$runID";
	}
    } elsif ( -e "$CMRUN/$runID" ) {            # exists and is not a link
	die "Can't create a link, a file/dir $CMRUN/$runID is on the way";
    } else {
        symlink "$RunDir", "$CMRUN/$runID" or die "Can't create link from $RunDir to $CMRUN/$runID";
    }
}

## Make sure that we really can write to dir/link $runID
if ( ! ( -d $runID && -r $runID && -w $runID && -x $runID ) ) {
    print "Couldn't create the link $runID\n";
    print "or run directory has wrong permissions. Aborting setup.\n";
    exit 1;
}

## If NetCDF template file is present copy it to run directory
if ( -s $netcdf_template_file ) {
    `cp $netcdf_template_file $runID/nctemp`;
}

## Switching to Run directory
chdir "$runID" or die "Can't chdir to $runID : $!\n";

open PRT, ">$runID".".PRT" or die "can't open ${runID}.PRT for writing\n";
print PRT "0Run $runID\n";

## Copy executable from "$runID"_bin directory to $runID
if ( ! -f "$DeckDir/$runID"."_bin/$runID.exe" ) {
    print "$DeckDir/$runID"."_bin/$runID.exe","not found\n";
    exit 1;
}

`cp -f "$DeckDir/$runID"_bin/$runID.exe $runID.exe`;
chmod 0777 & $umask_inv, "$runID.exe";

open RFILE, "$DeckDir/$rfile" or die "can't open $DeckDir/$rfile";

## Check signature
$_ = <RFILE>;
if ( ! /^\s*$runID/ ) {
    #print "inconsistent Naming: $rfile is not $_\n" ;
    #exit 1;
    #if first line doesn't start from runID add (prepend) it
    $_ = $runID." ".$_;
}
print PRT " $_";

## Read rfile until "Data input files:" is encountered
$run_options = 0;
while (<RFILE>) {
    print PRT " $_";
    if ( /^\s*Run *Options/ ) { $run_options = 1; }
    if ( /^\s*STACKSIZE *= *(\d+)/ && $run_options ) { $MIN_STACK = $1; }
    if ( /Data *input *files/ ) { last; }
}

print "Requested Stack  = $MIN_STACK\n" if ( $MIN_STACK ) ;

## Use the subsequent information to create the link and unlink files
while (<RFILE>) {
    if ( /Label *and *Namelist/ ) { last; }
    print PRT " $_";
    s/!.*//;
    push @data_files, /([\w.+_-]+\s*=\s*[\w.\/+_-]+)/g;
}

open RUNIDLN, ">$runID"."ln" or die "can't open ${runID}ln for writing\n";
open RUNIDULN, ">$runID"."uln" or die "can't open ${runID}ln for writing\n";

$flag_missing_data_files = 0;

foreach $_ ( @data_files ) {
    my $full_dest = "";
    ($name, $dest) = split /\s*=\s*/;
    if ( $dest !~ /^\// ) { 
	my $dir="";
	foreach $dir (split /:/, $GCMSEARCHPATH) { #/ - emacs fix
	    if ( -e "$dir/$dest" ) {
		$full_dest = "$dir/$dest";
		last;
	    }
	}
    } else {
	$full_dest = $dest;
    }
    if ( ! $full_dest || ! -e "$full_dest" ) {
	print "$dest not found\n";
	$flag_missing_data_files = 1;
	next;
    }
    print RUNIDLN "ln -fsn $full_dest $name\n";
    print RUNIDULN "rm $name\n";
    die "\n*** duplicate definition for $name file - please fix!\n\n" 
	if $data_file_hash{$name};
    $data_file_hash{$name} = $full_dest;
    push @data_file_names, $name;
}

if ( $flag_missing_data_files ) { exit 1; }

close RUNIDLN;
close RUNIDULN;

open RUNTIMEOPTS, ">runtime_opts" or die "can't create runtime_opts\n";
## Architecture-dependent settings
## Linux setting is only for Intel 7 compilers, but harmless otherwise
$uname = `uname`;
if ( $uname =~ /IRIX64/ ) {
    print RUNTIMEOPTS "export PAGESIZE_DATA=64 PAGESIZE_STACK=64\n";
} elsif ( $uname =~ /AIX/ ) {
    print RUNTIMEOPTS "export XLFRTEOPTS=NAMELIST=OLD\n";
} elsif ( $uname =~ /Linux/ ) {
    print RUNTIMEOPTS "export F_UFMTENDIAN=big\n";
    print RUNTIMEOPTS "export G95_ENDIAN=BIG\n";   #needed for g95
}
print RUNTIMEOPTS <<EOF;
    if [ `ulimit -s` != 'unlimited' ] ; then
      if [ `ulimit -s` -lt $MIN_STACK ] ; then
        ulimit -s $MIN_STACK || \\
          echo "!!! Could not set required stack size !!!" ;
        echo "current stack: `ulimit -s`"
      fi
    fi
EOF
close RUNTIMEOPTS;

## Architecture-dependent commands to start MPI runs
# default commands
$mpi_start = "";
$mpi_run = $MPIBIN."mpirun \$MPI_FLAGS -np \$NP ";
$mpi_stop = "";

if ( $MPIDISTR =~ /mvapich2/ ) {
    $mpi_run = $MPIBIN."mpirun_rsh -np \$NP -hostfile \$PBS_NODEFILE ";
} elsif ( $MPIDISTR =~ /SCALI/ ) {
    $mpi_run = $MPIBIN."mpirun \$MPI_FLAGS -np \$NP -inherit_limits ";
} elsif ( $MPIDISTR =~ /openmpi/ ) {
    if ( -f $MPIBIN."openmpirun" ) { # hack for Macports OpenMPI
	$mpi_run = $MPIBIN."openmpirun \$MPI_FLAGS -np \$NP";
    } else {
	$mpi_run = $MPIBIN."mpirun \$MPI_FLAGS -np \$NP"; # --mca btl_openib_warn_no_hca_params_found 0 ";
    }
}

# special cases for old architectures
if ( $uname =~ /OSF1/ && ! $MPIDISTR) {
    $mpi_run = "prun -s -n \$NP ";
}

# special cases for particulat machines
if ( $LOCATION =~ /Pleiades/i ) {
    if ( $MPIDISTR =~ /intel/ ) {
	$mpi_start = "mpdboot --file=\$PBS_NODEFILE --ncpus=1 --totalnum=`cat \$PBS_NODEFILE  | sort -u | wc -l` --ifhn=`head -1 \$PBS_NODEFILE` --rsh=ssh --mpd=`which mpd` --ordered";
	$mpi_run = "mpiexec \$MPI_FLAGS -np \$NP";
	$mpi_stop = "mpdallexit";
    } else {
        $mpi_run = "mpiexec \$MPI_FLAGS -np \$NP";
    }
}

# overwrite the command with the one from modelErc if present
if ( $MPIRUN_COMMAND ) {
    $mpi_run = $MPIRUN_COMMAND." ";
}

chmod 0777 & $umask_inv, "${runID}ln", "${runID}uln", "runtime_opts";

open I, ">I" or die "can't open 'I' for writing\n";
## Check signature
$_ = <RFILE>;
if ( ! /^\s*$runID[ (]/ ) {
    #print "inconsistent Naming: $rfile is not start of $_\n" ;
    #exit 1;
    $_ = $runID." ".$_;
}
print I;
$_ = <RFILE>;
print I;

## Hack: if in old format try to fix it
$istr = "";
$str_cold = "";
while (<RFILE>) {
    chop;
    s/!.*//g;
    s/^\s*/ /;
    s/\s*$//;
    next if ( ! $_ ); 
    #print I "$_\n";
    if ( /^\s*&&END_PARAMETERS/ ) {
	#foreach $key (sort (keys %data_file_hash)) {
	foreach $key (@data_file_names) {
	    $istr .= " _file_$key='$data_file_hash{$key}'\n";
	}
    }
    if ( /^\s*ISTART/ ) {
	$str_cold .= "$_\n";
    } else {
	$istr .= "$_\n";
    }
}

print I "$istr\n";

if ( $str_cold ) {
    print I " &INPUTZ_cold\n";
    print I "$str_cold";
    print I " /\n";
}

close I;
close PRT;
chmod 0666 & $umask_inv, "I", "${runID}.PRT";

if ( $mpi ) {
    $run_command = $mpi_run;
} else {
    $run_command = "";
}


## Create executable script file RUNID
open RUNID, ">$runID" or die "can't open $runID for writing\n";
print RUNID <<EOF;
\#!/bin/sh
    trap '' TERM
    PRTFILE=${runID}.PRT
    IFILE="I"
    NP="\$MP_SET_NUMTHREADS"
    debug=0
    opts=
    touch_ifile=0
    if [ "\$NP"x = x ] ; then NP=1; fi
    if [ "\$DEBUG_COMMAND"x = x ] ; then DEBUG_COMMAND="xterm -e gdb --args"; fi
    while [ \$\# -ge 1 ] ; do
      OPT=\$1 ; shift
      case \$OPT in
        -q)
            PRTFILE='/dev/null'
            ;;
        -l)
            PRTFILE="\$1" ; shift
            ;;
        -i)
            IFILE="\$1" ; shift
            ;;
        -np)
            NP="\$1" ; shift
            ;;
        -d)
            debug=1
            ;;
        -cold-restart)
            opts="\$opts -cold-restart" ; touch_ifile=1
            ;;
        --time)
            opts="\$opts --time \$1" ; shift
            ;;
         *)
            echo "Warning: wrong option ignored: \$OPT"
            ;;
      esac
    done
    umask $umask_str
    rc=99
    if [ -s run_status ] ; then rc=`head -1 run_status` ; fi
    if [ \$rc -eq 13 ] ; then
      if [ `find run_status -newer I` ] ; then
        echo 'run seems to have finished'
        exit 1; fi
    fi
    if [ -f lock ] ; then
      echo 'lock file present - aborting' ; exit 1 ; fi
    touch lock
    echo '-99' > run_status
    echo 'INPUT not yet completed' >> run_status
    . ./runtime_opts
    if [ -s nctemp ] ; then $NETCDFBIN/ncgen -b -o nctemp.nc nctemp ; fi
    ./${runID}ln
    $mpi_start
    if [ \$debug -eq 1 ] ; then
      $run_command \$DEBUG_COMMAND ./${runID}.exe -i ./\$IFILE \$opts
    else
      $run_command ./${runID}.exe -i ./\$IFILE \$opts > \$PRTFILE
    fi
    $mpi_stop
    rc=`head -1 run_status`
    ./${runID}uln
    rm -f lock
    if [ \$touch_ifile -eq 1 ] ; then sleep 1 ; touch \$IFILE ; fi
    exit \$rc
EOF
close RUNID;
chmod 0777 & $umask_inv, $runID;
## end of RUNID script
`ln -sf $runID E`;

## save the list of current modules
`echo \$LOADEDMODULES | sed -s 's/:/ /g;' > modules`;

## setup finished normally
exit 0 ;

