editRundeck()
# -------------------------------------------------------------------
{
# Adapted from R. Ruedy's script add_regr_test_params
# This function enables regression testing of restarts. It modifies
# the RUNDECK by adding an NDISK line right before the &&END_PARAMETERS
# and an end time near the end of the RUNDECK.

  local deck=$1".R"
  local ndisk=$2
  local datee=$3
  local houre=$4

  echo $deck
  ndisk_line="ndisk=$ndisk"
# If there is no fifth argument (default) then:
  if [ -z "$5" ]; then
    end_hour_line=" DATEE=$datee, HOURE=$houre,"
# else modify MONTHE (and potentially YEARE):
  else
    # exclude commented (!) lines, get first match, print year value
    year=`grep  "^[^\!]" $deck | grep -m1 YEARE | awk -F= '{print $2}' | awk -F, '{print $1}'`
    # exclude commented (!) lines, get first match, print month value
    month=`grep "^[^\!]" $deck | grep -m1 YEARE | awk -F= '{print $3}' | awk -F, '{print $1}'`
    if [ $month -ge 11 ]; then
      let "year=year+1"
      if [ $month -eq 12 ]; then
        let "month=2"
      else
        let "month=1"
      fi
    else
      let "month=month+2"
    fi
    end_hour_line=" YEARE=$year, MONTHE=$month, DATEE=$datee, HOURE=$houre,"
  fi
  eof1='&&END_PARAMETERS'
  eof2='/'

  a=$( grep -n ${eof1}     $deck | head -1 ) ; n1=${a%%:*}
  a=$( grep -n ^${eof2} $deck | head -1 ) ; n2=${a%%:*}

  cp ${deck} templ
  head -$(( n1-1 )) templ                   > ${deck}
  echo "${ndisk_line}"                      >> ${deck}
  tail +${n1} templ | head -$(( n2 - n1 ))  >> ${deck}
  echo "${end_hour_line}"                   >> ${deck}
  if [[ ${deck} =~ "ENINT" ]]; then
     echo " ISTART=8, ${end_hour_line}"         >> ${deck}
  else
     echo " ISTART=2, ${end_hour_line}"         >> ${deck}
  fi
  tail +${n2} templ                         >> ${deck}
}

editRundeck $1 $2 $3 $4 $5
