#! /bin/sh
# ward_driver.sh --
########################################################################

omega="$1"
shift

models="qed qcd sym sm sm_top_anom"

modules=""

while read module threshold n roots model unphysical mode process; do 

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
      eval n_$module=$n
      eval roots_$module=$roots
      eval process_$module="'$process'"
      ########################################################################

      ########################################################################

      omega_bin="`echo $omega | sed s/%%%/$model/g`"
      # echo "running $omega_bin -$mode '$process'" 1>&2
      $omega_bin "$@" \
        -target:parameter_module parameters_$model \
        -target:module amplitude_ward_physical_$module \
        -$mode "$process" 2>/dev/null
      $omega_bin "$@" \
        -target:parameter_module parameters_$model \
        -target:module amplitude_ward_unphysical_$module \
        -$mode "$process" -unphysical $unphysical 2>/dev/null
    ;;
  esac

done
########################################################################

for module in $modules; do

    for mode in physical unphysical; do

cat <<EOF
module interface_ward_${mode}_${module}
  use omega_interface
  use amplitude_ward_${mode}_${module}
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
    p%external_masses => external_masses
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
end module interface_ward_${mode}_${module}

EOF

    done

done

########################################################################

cat <<EOF
program ward_driver
  use kinds
  use ward_lib
EOF

for module in $modules; do
    for mode in physical unphysical; do
cat <<EOF
  use interface_ward_${mode}_${module}, load_${mode}_${module} => load
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
  real(kind=default), parameter :: THRESHOLD = 0.8
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
eval roots="\${roots_$module}"

cat <<EOF
  print *, "checking process '$process'"
  call check (load_physical_$module (), load_unphysical_$module (), &
              roots = real ($roots, kind=default), &
              threshold = real ($threshold, kind=default), &
              n = $n, seed = SEED, &
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
end program ward_driver
EOF

exit 0
