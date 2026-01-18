#! /bin/sh
# ects_driver.sh --
########################################################################

omega="$1"
shift

modules=""
modules_count=0

########################################################################
for ects_name in "$@"; do

 case $ects_name in

  # Allow ourselves as an argument and skip it.
  # This simplifies automake rules ...
  *.sh) :
    ;;

  *)
  ########################################################################
  module="`basename $ects_name .ects`"
  modules="$modules $module"
  model="`sed -n '/^!ects:model */s///p' $ects_name`"
  process="`sed -n '/^!ects:process */s///p' $ects_name`"
  if grep '^!ects:scatter' $ects_name >/dev/null 2>&1; then
    mode=scatter
  elif grep '^!ects:decay' $ects_name >/dev/null 2>&1; then
    mode=decay
  else
    mode=scatter
  fi
  ########################################################################

  omega_bin="`echo $omega | sed s/%%%/$model/g`"
  # echo "running $omega_bin -$mode '$process'" 1>&2
  $omega_bin \
    -target:parameter_module parameters_$model \
    -target:module amplitude_color_flows_$module \
    $omega_opts -$mode "$process" 2>/dev/null

  cat <<EOF
module expected_color_flows_$module
  use kinds
  use omega_color, only: OCF => omega_color_factor
  implicit none
  private
  integer :: i
EOF

  cat $ects_name

  cat <<EOF
  public :: cflows, gflags, cfactors
end module expected_color_flows_$module
EOF

  cat <<EOF
module test_color_flows_$module
  use color_test_lib
  use amplitude_color_flows_$module
  use expected_color_flows_$module
  implicit none
contains
  subroutine test (ok)
    logical, intent(out) :: ok
    call compare_color_flows ("$process", cflows, gflags, cfactors, &
         number_particles_in, number_particles_out, &
         number_color_indices, number_color_flows, color_flows, &
         number_color_factors, color_factors, ok)
  end subroutine test
end module test_color_flows_$module
EOF

  ;;
 esac

done
########################################################################


########################################################################
cat <<EOF
program main_$module
EOF

for module in $modules; do
  echo "  use test_color_flows_$module, test_$module => test"
done

cat <<EOF
  implicit none
  logical :: ok
  integer :: pass, fail
  pass = 0
  fail = 0
EOF

for module in $modules; do
  cat <<EOF
  call test_$module (ok)
  if (ok) then
     pass = pass + 1
  else
     fail = fail + 1
  end if
EOF
done

cat <<EOF
  print *, pass, "tests passed, ", fail, " tests failed."
  if (fail .gt. 0) then
    stop 1
  else 
    stop 0
  end if
end program main_$module
EOF

exit 0
