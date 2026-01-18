#! /bin/sh
# fermi_driver.sh --
########################################################################

omega="$1"
shift

models="qed qcd sym sm sm_top_anom mssm"

modules=""

########################################################################
while read prefix threshold abs_threshold n roots model i j eps mode process; do 

  case $prefix in

   '#'*) # skip comments
     ;;

   '')   # skip empty lines
     ;;

   '!'*) break
     ;;

    *)
      ########################################################################
      module=${prefix}_${i}_${j}
      modules="$modules $module"
      eval threshold_$module=$threshold
      eval abs_threshold_$module=$abs_threshold
      eval n_$module=$n
      eval i_$module=$i
      eval j_$module=$j
      eval eps_$module=$eps
      eval roots_$module=$roots
      eval process_$module="'$process'"
      ########################################################################


      omega_bin="`echo $omega | sed s/%%%/$model/g`"
      # echo "running $omega_bin -$mode '$process'" 1>&2
      $omega_bin "$@" \
        -target:parameter_module parameters_$model \
        -target:module amplitude_fermi_$module \
        -$mode "$process" 2>/dev/null
    ;;
  esac

done
########################################################################

for module in $modules; do

cat <<EOF
module interface_fermi_${module}
  use omega_interface
  use amplitude_fermi_${module}
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
end module interface_fermi_${module}

EOF

done

########################################################################

cat <<EOF
program fermi_driver
  use kinds
  use fermi_lib
EOF

for module in $modules; do
cat <<EOF
  use interface_fermi_${module}, load_${module} => load
EOF
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
  real(kind=default), parameter :: ABS_THRESHOLD = 1.0E-11
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
eval i="\${i_$module}"
eval j="\${j_$module}"
eval eps="\${eps_$module}"
eval threshold="\${threshold_$module}"
eval abs_threshold="\${abs_threshold_$module}"
eval roots="\${roots_$module}"

cat <<EOF
  print *, "checking process '$process' ($i <=> $j)"
  call check (load_$module (), i = $i, j = $j, eps = $eps, &
              roots = real ($roots, kind=default), &
              threshold = real ($threshold, kind=default), &
              abs_threshold = real ($abs_threshold, kind=default), &
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
end program fermi_driver
EOF

exit 0
