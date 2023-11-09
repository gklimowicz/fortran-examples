#!/bin/ksh 
SCALE=~mkelley5/bin/scaleacc_112511

# expects full years of CACHED_SUBDD subdd output files in DATA/ subdir

usage () {
  echo
  echo "Usage: $0 [-L -n -S] runName year1 year2 [yearOffset]"
  echo
  echo "  This script converts non-instantaneous 3D MRO3 and 2D p_surf daily CACHED_SUBDD"
  echo "subdd output files into daily NINT model ozone input for only the model levels"
  echo "of the subdd files. (I.e. pasting above the GCM top and any unit conversions"
  echo "are done online in the NINT run that will read these files.) This script expects"
  echo "your subdd files to be in the ./DATA subdirectory and (for now) that there are"
  echo "full years present. By default, one file per year is written. These are also"
  echo "linked to YYYY.nc filenames if the -L option is used. However, if the -S flag"
  echo "is employed a single file is created instead."
  echo "  If the -L or the -S options are chosen then a 4th argument for year offset is"
  echo "required (for example, 150 if run year 1851 is to represent real year 2001.)"
  echo " "
  echo "-L option to symbolicly link resulting files to YYYY.nc files for model reading"
  echo "   after offsetting the year by an integer from argument 4"
  echo " "
  echo "-S make a single model input file instead of one per year. In this case the 4th"
  echo "   argument offset is used for the time dimension meta data instead of linking."
  echo " "
  echo "-n no conversion (e.g. it was done alredy). Only do the linking (if -L is on)."
  echo
  return
}

# functions first:

scaleAndExtract() {
  file=${year}${MNUM[$m]}${DNUM[$d]}.subdd${run}.nc
  ofile=${year}${MNUM[$m]}${DNUM[$d]}.taijlh48${run}.nc
  pfile=${year}${MNUM[$m]}${DNUM[$d]}.aijh48${run}.nc
  $SCALE DATA/${file} all > /dev/null
  # After output file is scaled, extract the two variables needed and
  # concatenate into one file per year:
  for v in MRO3 p_surf Ozone ; do
    if [[ -e ${v}_now.nc ]] ; then rm ${v}_now.nc ; fi
  done
  ncks -v MRO3 ${ofile} MRO3_now.nc ; rm $ofile
  ncks -v p_surf ${pfile} p_surf_now.nc ; rm $pfile
  ncks p_surf_now.nc Ozone_now.nc ; rm p_surf_now.nc
  ncks -A MRO3_now.nc Ozone_now.nc ; rm MRO3_now.nc
  return
}

appendToFile() {
  thisFile=$1
  if [[ $fileIsNew == "yes" ]] ; then
    mv Ozone_now.nc ${thisFile}
    fileIsNew="no"
  else
    mv ${thisFile} temp_${thisFile}
    # append to output file:
    ncrcat -h temp_${thisFile} Ozone_now.nc ${thisFile} ; rm temp_${thisFile} 
  fi
  return
}

# execution begins here:

# set defaults and read option flags:
singleFile="no"
doLinks="no"
doConvert="yes"
while getopts "LnS" option;
do
 case $option in
  L) ; doLinks="yes" ;;
  n) ; doConvert="no" ;;
  S) ; singleFile="yes" ;;
  *) ; usage ; echo "invalid option -$option" ; exit ;;
 esac
done
shift $(($OPTIND - 1))

if [[ $# -lt 3 || ( $doLinks == "yes" && $# -lt 4 ) || ( $singleFile == "yes" && $# -lt 4 ) ]] ; then
  usage
  echo "INSUFFICIENT ARGUMENTS."
  exit 13
fi

if [[ $doConvert == "yes" ]] ; then
  echo "Host appears to be: $HOST."
  echo "Important! Should be run on backend (qsubI) node because of himem version of scaleacc!"
  echo "Do control-C to quit; return to continue."
  read nothing
fi

if [[ ! -d OUT ]] ; then mkdir OUT ; fi

# assume no leap years:
set -A MNUM 01 02 03 04 05 06 07 08 09 10 11 12
set -A MDAY 31 28 31 30 31 30 31 31 30 31 30 31
set -A DNUM 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

# read the arguments:
run=$1
y1=$2
y2=$3
if [[ $doLinks == "yes" || $singleFile == "yes" ]] ; then yO=$4 ; fi # offset year by yO

if [[ $singleFile == "yes" ]] ; then
  ### PUT ALL DATA INTO SINGLE MODEL INPUT FILE (timestream method 1) ###

  daycount=0
  daycountM1=-1
  year=$y1
  offsettedYear1=`echo "$y1 + $yO" | bc`
  offsettedYear2=`echo "$y2 + $yO" | bc`
  bigFile=Daily_Ozone_from_${run}_years_${offsettedYear1}-${offsettedYear2}.nc
  if [[ -e ${bigFile} ]] ; then rm ${bigFile} ; fi
  echo "Doing: ${bigFile}"
  fileIsNew="yes"
 
  while [[ $year -le $y2 ]]; do # loop over years
    offsettedYear=`echo "$year + $yO" | bc`
    m=0 
    while [[ $m -le 11 ]] ; do # loop over months (-1 because indicies)
      echo "  Month=${MNUM[$m]} Year=${year}"
      d=0
      while [[ $d -le $( echo "${MDAY[$m]} - 1" | bc ) ]] ; do # loop over days (-1 because indicies)
        let daycount+=1
        let daycountM1+=1
        scaleAndExtract
        appendToFile ${bigFile}
        # append to time dimension:
        tuple_edit $bigFile $daycount $daycount time time ${daycountM1}.5 > /dev/null 
        let d+=1
      done # days
      let m+=1
    done # months
    let year+=1
  done # years
  # change the units of the time coordinate variable:
  ncatted -O -a units,time,o,c,"days since ${offsettedYear1}-01-01" ${bigFile}
  mv -f ${bigFile} OUT/
  rm -f Ozone_now.nc

else 
  ### PUT DATA INTO ONE FILE PER YEAR (timestream method 2) ###

  year=$y1
  while [[ $year -le $y2 ]]; do # loop years
    yearFile=Daily_Ozone_from_${run}_year_${year}.nc
    if [[ $doConvert == "yes" ]] ; then
      fileIsNew="yes"
      if [[ -e ${yearFile} ]] ; then rm ${yearFile} ; fi
      echo "Doing: ${yearFile}"
      m=0
      while [[ $m -le 11 ]] ; do # loop months (-1 because indicies)
        echo "  Month=${MNUM[$m]} Year=${year}"
        d=0 
        while [[ $d -le $( echo "${MDAY[$m]} - 1" | bc ) ]] ; do # loop days (-1 because indicies)
          scaleAndExtract
          appendToFile ${yearFile}
          let d+=1
        done # days
        let m+=1
      done # months
      mv -f ${yearFile} OUT/
    fi
    if [[ $doLinks == "yes" ]] ; then
      offsettedYear=`echo "$year + $yO" | bc`
      ln -s OUT/${yearFile} ./${offsettedYear}.nc
    fi
    let year+=1
  done # years
  if [[ $doConvert == "yes" ]]; then rm Ozone_now.nc ; fi

fi # single/multiple files choice

