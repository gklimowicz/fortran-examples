#! /bin/sh
# compare_driver.sh --
########################################################################

tag=VM
odefault="$1"
ovm="$2"
oparams="$2 $3"
shift 3

modules=""

########################################################################
while read module threshold abs_threshold n roots model mode process; do

  case $module in

   '#'*) # skip comments
     ;;

   '')   # skip empty lines
     ;;

   '!'*) break
     ;;

    *)
      modules="$modules $module"
      eval threshold_$module=$threshold
      eval abs_threshold_$module=$abs_threshold
      eval n_$module=$n
      eval roots_$module=$roots
      eval process_$module="'$process'"
      eval model_$module="'$model'"
    ;;
  esac

done
########################################################################

for module in $modules; do

  eval bc_file="\${bc_file_$module}"
  eval model="\${model_$module}"
  ovp="`echo $oparams | sed s/%%%/$model/g`"
  ovp="`echo $ovp | sed s/%%/amplitude_compare_${tag}_v2_${module}/g`"
  ovp="`echo $ovp | sed s/%/${module}.hbc/g`"
  #echo "writing wrapper file $params_file with '$ovp'" 1>&2
  $ovp 2>/dev/null

done

########################################################################

exit 0
