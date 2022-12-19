#! /bin/sh
# compare_driver.sh --
########################################################################

tag="$1"
omega_v1="$2"
omega_v2="$3"
shift 3

models="qed qcd sym"

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
      ########################################################################
      modules="$modules $module"
      eval threshold_$module=$threshold
      eval abs_threshold_$module=$abs_threshold
      eval n_$module=$n
      eval roots_$module=$roots
      eval process_$module="'$process'"
      ########################################################################

      ov1="`echo $omega_v1 | sed s/%%%/$model/g`"
      # echo "running $ov1 -$mode '$process'" 1>&2
      $ov1  "$@" \
        -target:parameter_module parameters_$model \
        -target:module amplitude_compare_${tag}_v1_${module} \
        -$mode "$process" 2>/dev/null
      ov2="`echo $omega_v2 | sed s/%%%/$model/g`"
      # echo "running $ov2 -$mode '$process'" 1>&2
      $ov2 "$@" \
        -target:parameter_module parameters_$model \
        -target:module amplitude_compare_${tag}_v2_${module} \
        -$mode "$process" 2>/dev/null
    ;;
  esac

done
########################################################################

for module in $modules; do

    for mode in v1 v2; do

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
program compare_driver
  use kinds
  use compare_lib
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
  integer, parameter :: N = 1000
  real(kind=default), parameter :: ROOTS = 1000
  integer, parameter :: SEED = 42
  integer :: failures, attempts, failed_processes, attempted_processes
  failed_processes = 0
  attempted_processes = 0
EOF

for model in $models; do
cat <<EOF
  call init_parameters_$model ()
EOF
done

for module in $modules; do

eval process="\${process_$module}"
eval n="\${n_$module}"
eval threshold="\${threshold_$module}"
eval abs_threshold="\${abs_threshold_$module}"
eval roots="\${roots_$module}"

cat <<EOF
  print *, "checking process '$process'"
  call check (load_v1_$module (), load_v2_$module (), &
              roots = real ($roots, kind=default), &
              threshold = real ($threshold, kind=default), &
              n = $n, seed = SEED, &
              abs_threshold = real ($abs_threshold, kind=default), &
              failures = failures, attempts = attempts)
  if (failures > 0) then
     print *, failures, " failures in ", attempts, " attempts"
     failed_processes = failed_processes + 1
  end if
EOF
done

cat <<EOF
  if (failed_processes > 0) then
     print *, failed_processes, " failed processes in ", attempted_processes, " attempts"
     stop 1
  end if
end program compare_driver
EOF

exit 0
