# Usage:
#
# extract_scm.sh RUN_NAME.R
#
# For SCM rundeck RUN_NAME.R specifying geographic location lon_targ, lat_targ
# in its parameters section, this script extracts data for that location
# from the set of rundeck input files having the pattern
#
# SHORTNAME=/.+/extractions/filename.nc
# or
# SHORTNAME=/.+/extractions/some/directory/name/YYYY.nc
#
# See contextualized usage notes in the input files section of templates/SCM.R
#
# This script is very rough and will be rewritten in either perl or python.
# NCO is currently used for the actual extraction, but a special-purpose
# utility might be written if it proves more convenient w.r.t.
# issues of spatial interpolation and arbitrary input grids.
#



rundeck=$1

rcfile=${HOME}/.modelErc
if [[ ! -s ${rcfile} ]]; then
  echo "cannot find ${rcfile}"
  exit
fi

#
# Obtain the set of default input paths
#
inpdirs=`sed 's/\#.*//' ${rcfile} | grep 'GCMSEARCHPATH=' | sed 's/=/ /' | awk '{print $2}' | sed 's/:/ /g'`

#
# Obtain the desired lon/lat information from the rundeck.
#

lon_targ=`grep 'SCM_lon=' $rundeck | sed 's/!.*//' | tail -1 | sed 's/=/ /' | awk '{print $2}'`
lat_targ=`grep 'SCM_lat=' $rundeck | sed 's/!.*//' | tail -1 | sed 's/=/ /' | awk '{print $2}'`

if [[ $lon_targ == '' ]]; then echo "missing lon_targ"; exit; fi
badlon=`echo $lon_targ | awk '$1<-180.||$1>180.{print "badlon"}'`
if [[ $badlon == "badlon" ]]; then echo "bad lon_targ"; exit; fi

if [[ $lat_targ == '' ]]; then echo "missing lat_targ"; exit; fi
badlat=`echo $lat_targ | awk '$1<-90.||$1>90.{print "badlat"}'`
if [[ $badlat == "badlat" ]]; then echo "bad lat_targ"; exit; fi

#
# Before extracting any data, compose the list of full input
# paths and perform sanity checks.
#
flist0=`egrep '=/.+/extractions/' $rundeck | sed 's/!.*//' | sed 's/=/ /' | awk '{print $2}'`

iolist=''
firstbase=''

for fname in $flist0; do

parts=`echo $fname | sed 's/\/extractions\//\/extractions /'`
outbase=`echo $parts | awk '{print $1}'`
bname=`echo $parts | awk '{print $2}'`

inpfile=''
for inpdir in $inpdirs
do
  inpfilex=${inpdir}"/"${bname}
  if [[ -s ${inpfilex} ]]; then
    inpfile=${inpfilex}
    break
  fi
done
if [[ $inpfile == '' ]]; then
  echo "parent file for ${fname} does not exist in GCMSEARCHPATH: exiting"
  exit
fi

if [[ $firstbase == '' ]]; then firstbase=$outbase; fi
if [[ $outbase != $firstbase ]]; then
  echo "Please change rundeck filename:"
  echo $fname
  echo "to:"
  echo "${firstbase}/${bname}"
  echo "or vice versa, so that all extracted files for the same"
  echo "lon_targ,lat_targ are placed into the same directory.  Exiting."
  exit
fi

if [[ ${inpfile} == ${fname} ]]; then
  echo "input and output file paths are identical for:"
  echo ${fname}
  echo "Exiting."
  exit
fi

iolist=${iolist}" ${inpfile}::::${fname}"

done

#
# Eliminate duplicate files from the list
#
iolist=`echo $iolist | tr ' ' '\n' | sort | uniq`


#
# Expand directory names containing files of the form YYYY.nc
#
iolist0=${iolist}
iolist=''
for iopair in $iolist0; do

inpfile=`echo $iopair | sed 's/::::/ /' | awk '{print $1}'`
outfile=`echo $iopair | sed 's/::::/ /' | awk '{print $2}'`

if [[ ${inpfile} == ${inpfile%.nc} ]]; then
#if [[ -d ${inpfile} ]]; then
  for yrfile in `ls ${inpfile}/[0-9][0-9][0-9][0-9].nc`; do
    yr=`basename ${yrfile}`
    iolist=${iolist}" ${yrfile}::::${outfile}/${yr}"
  done
else
  iolist=${iolist}" "${iopair}
fi

done

#
# Check if netcdf files have the requisite coordinate info
#
for iopair in $iolist; do

inpfile=`echo $iopair | sed 's/::::/ /' | awk '{print $1}'`
outfile=`echo $iopair | sed 's/::::/ /' | awk '{print $2}'`

lonandlat=`ncdump -h $inpfile | egrep " lon\(lon\)| lat\(lat\)" | wc -l`
if [[ $lonandlat -ne 2 ]]; then
  echo "file $inpfile"
  echo "does not contain lon and lat coordinate variables.  Exiting."
  exit
fi

# restart files are not dimensioned by lon and lat.  Look for requisite auxiliary info.
iandj=`ncdump -h $inpfile | egrep "i_of_lon\(lon\)|j_of_lat\(lat\)" | wc -l`
if [[ $iandj -ne 0 ]]; then
if [[ $iandj -ne 2 ]]; then
  echo "restart file $inpfile"
  echo "does not contain i,j index info.  Exiting."
  exit
fi
fi

done

#
# Perform the extraction
#

for iopair in $iolist; do

inpfile=`echo $iopair | sed 's/::::/ /' | awk '{print $1}'`
outfile=`echo $iopair | sed 's/::::/ /' | awk '{print $2}'`
dname=`dirname ${outfile}`

if [[ ! -s $dname ]]; then
  echo "${dname} does not exist: creating"
  mkdir -p $dname
fi

if [[ -s $dname ]]; then

cmd="ncks -O -a -d lon,${lon_targ} -d lat,${lat_targ} $inpfile $outfile"
echo $cmd
$cmd

# restart files are not dimensioned by lon and lat.  extract by i,j indices
iandj=`ncdump -h $outfile | egrep "i_of_lon\(lon\)|j_of_lat\(lat\)" | wc -l`
if [[ $iandj -eq 2 ]]; then
  i_targ=`ncdump -v i_of_lon $outfile | grep "i_of_lon =" | awk '{print $3}'`
  j_targ=`ncdump -v j_of_lat $outfile | grep "j_of_lat =" | awk '{print $3}'`
  iname=`ncdump -h $outfile | egrep "lon:alt_dimname" | sed 's/"//g' | awk '{print $3}'`
  jname=`ncdump -h $outfile | egrep "lat:alt_dimname" | sed 's/"//g' | awk '{print $3}'`
  cmd="ncks -O -a -F -d ${iname},${i_targ} -d ${jname},${j_targ} $outfile $outfile"
  echo $cmd
  $cmd
fi

fi

done
