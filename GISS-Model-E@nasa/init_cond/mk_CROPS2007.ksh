#!/usr/bin/ksh

act=/usr/local/other/NCO/3.9.5_intel9.1.052/bin/ncks
if [[ ! -s $act ]]
then echo "this script uses <ncks> - netcdf's kitchen sink utility"
     echo "the current path $act seems wrong - please correct it"
     exit ; fi

if [[ $# -ne 2 ]]
then echo "Usage: $0 IMxJM topo_file"
     echo
     echo "the name of the resulting crops file will be CROPS2007_topo-Z"
     echo "          where topo-Z is topo_file without the leading Z"
     echo "examples: $0 72x46 Z72X46N.cor4_nocasp => CROPS2007_72X46N.cor4_nocasp"
     echo "          $0 144x90 Z144X90N_nocasp    => CROPS2007_144X90N_nocasp"
     exit ; fi

  ### input files (crops files are listed in $lst)
grid=$1 ; topo=$2 ; lst=crops_years_used.txt

  ### output file $out
topo_name=$( basename $topo ) ; out=CROPS2007_${topo_name#Z} ; # strip leading Z

  ### check command line arguments and existence of needed files
if [[ $grid != *[xX]* ]] ; then echo "$grid - not of the form IMxJM" ;    exit ; fi
if [[ ! -s $topo ]] ; then echo "topogr_file $topo not found" ;           exit ; fi
  echo
  echo "processing vegetation files archived on /archive/u/mpuma/Globalveg"
  echo "  the following needs to have been done before using this script:"
  echo "  copy them, unzip them and move them to the current directory"
  echo "checking whether all the needed files are present ..."
argnc='' ; argtxt=''
while read a
do if [[ ! -s $a ]]
   then echo ; echo "$a not found - no action" ; echo
        exit ; fi
   argnc="$argnc $a" ; argtxt="$argtxt ${a%nc}txt"
done < ${lst}

  echo "extracting the crops data as text files ..."
  echo
  echo "ignore the ncks messages - hit enter if they fill the screen and things stop"
  echo

for x in $argnc
do ${act} -d vegtype,11 $x | grep time | grep vegtype 2> /dev/null > ${x%nc}txt
done
  echo ; echo "done with extracting" ; echo

  echo "compiling regridding/rescaling program"
ifort -convert big_endian convert_crops.f -o convert_crops
  echo "converting ..."  ; convert_crops $grid $topo $argtxt
mv -f CROPS2007.ext $out        ; echo "created $out"
rm -f convert_crops $argtxt ; # clean up
