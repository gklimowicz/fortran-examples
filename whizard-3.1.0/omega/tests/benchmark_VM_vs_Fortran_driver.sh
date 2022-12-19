#! /bin/sh
# compare_driver.sh --
########################################################################

tag=VM
odefault="$1"
ovm="$2"
oparams="$2 $3"
shift 3

models="QED QCD SM"

modules=""

########################################################################
while read module n roots model mode process; do

  case $module in

   '#'*) # skip comments
     ;;

   '')   # skip empty lines
     ;;

   '!'*) break
     ;;

    *)
      ########################################################################
      modules="$modules $module"
      eval n_$module=$n
      eval roots_$module=$roots
      eval process_$module="'$process'"
      eval model_$module="'$model'"
      ########################################################################

      ov1="`echo $odefault | sed s/%%%/$model/g`"
      #echo "running $ov1 -$mode '$process'" 1>&2
      $ov1  "$@" \
        -target:parameter_module parameters_$model \
        -target:module amplitude_compare_${tag}_v1_${module} \
        -$mode "$process" 2>/dev/null
      ov2="`echo $ovm | sed s/%%%/$model/g`"
      #bc_file="`echo $process | sed 's/->/to/g' | sed 's/ /_/g'`".hbc
      bc_file=$module.hbc
      eval bc_file_$module=$bc_file
      #echo "running $ov2 -$mode '$process', saving to $bc_file" 1>&2
      $ov2 "$@" -$mode "$process" 2>/dev/null 1> $bc_file
    ;;
  esac

done
########################################################################

for module in $modules; do

    for mode in v1 v2; do

      if [ $mode == v2 ]; then
        init="call init()"
      else
        init=""
      fi

cat <<EOF
module interface_compare_${tag}_${mode}_${module}
  use omega_interface
  use amplitude_compare_${tag}_${mode}_${module}
  implicit none
  private
  public :: load
contains
  function load () result (p)
    type(omega_procedures) :: p
    $init
    p%number_particles_in => number_particles_in
    p%number_particles_out => number_particles_out
    p%number_spin_states => number_spin_states
    p%spin_states => spin_states
    p%number_flavor_states => number_flavor_states
    p%flavor_states => flavor_states
    p%number_color_indices => number_color_indices
    p%number_color_flows => number_color_flows
    p%color_flows => color_flows
    p%number_color_factors => number_color_factors
    p%color_factors => color_factors
    p%color_sum => color_sum
    p%new_event => new_event
    p%reset_helicity_selection => reset_helicity_selection
    p%is_allowed => is_allowed
    p%get_amplitude => get_amplitude
  end function load
end module interface_compare_${tag}_${mode}_${module}

EOF

    done

done

########################################################################

cat <<EOF
program benchmark_VM_vs_Fortran
  use kinds
  use benchmark_lib
EOF

for module in $modules; do
    for mode in v1 v2; do
cat <<EOF
  use interface_compare_${tag}_${mode}_${module}, load_${mode}_${module} => load
EOF
    done
done

for model in $models; do
cat <<EOF
  use parameters_$model, init_parameters_$model => init_parameters
EOF
done

cat <<EOF
  implicit none
  integer, parameter :: SEED = 42
  real(double) :: time1, time2
  write(*, "(A)") '#All times in ms/phasespacepoint. All colors and helicities computed.'
  write(*, "(A)") '#Process      Fortran                 VM         VM/Fortran'
EOF

for model in $models; do
cat <<EOF
  call init_parameters_$model ()
EOF
done

for module in $modules; do

eval process="\${process_$module}"
eval n="\${n_$module}"
eval roots="\${roots_$module}"

cat <<EOF
  call benchmark (time1, time2, load_v1_$module (), load_v2_$module (), &
              roots = real ($roots, kind=default), &
              n = $n, seed = SEED)
              write(*, "(A,F19.5,F19.5,F19.5)") "$module", time1, time2, time2 / time1
EOF
done

cat <<EOF
end program benchmark_VM_vs_Fortran
EOF

exit 0
