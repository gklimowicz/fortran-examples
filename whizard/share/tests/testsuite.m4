dnl $Id: testsuite.m4 5205 2014-01-30 14:17:39Z jr_reuter $
dnl ====================================================================
dnl
define(`TS__PRC_NUMBER',`0')dnl
define(`TS__DIV_COMPILE',`1')dnl
define(`TS__DIV_INTEGRATE',`2')dnl
define(`TS__DIV_COMPARE',`3')dnl
dnl
define(`TS__PRC_PREFIX',`prc')dnl
dnl
dnl ====================================================================
dnl
dnl PROCESS (IN-PARTICLES, OUT-PARTICLES, SQRTS, SIGMA, ACCURACY, TOLERANCE)
dnl
dnl --------------------------------------------------------------------
define(`PROCESS',`dnl
dnl ====================================================================
dnl
dnl Construct a unique process name from a counter
dnl 
define(`TS__PRC_NUMBER',incr(TS__PRC_NUMBER))dnl
define(`TS__PRC_NAME',TS__PRC_PREFIX`'TS__PRC_NUMBER)dnl
dnl 
dnl Store the arguments in descriptive macros (we need the
dnl versions with stripped double quotes for use in printf).
dnl 
define(`TS__PRC_IN_PARTICLES',`$1')dnl
define(`TS__PRC_IN_PARTICLES_UNQUOTED',`translit(`$1',`"')')dnl
define(`TS__PRC_OUT_PARTICLES',`$2')dnl
define(`TS__PRC_OUT_PARTICLES_UNQUOTED',`translit(`$2',`"')')dnl
define(`TS__PRC_SQRTS',`$3')dnl
define(`TS__PRC_INTEGRAL',`$4')dnl
define(`TS__PRC_ACCURACY',`$5')dnl
define(`TS__PRC_TOLERANCE',`$6')dnl
define(`TS__PRC_DESCRIPTION',`TS__PRC_IN_PARTICLES_UNQUOTED => TS__PRC_OUT_PARTICLES_UNQUOTED @ sqrt(s) = TS__PRC_SQRTS')dnl
dnl
dnl Process definition
dnl
divert(TS__DIV_COMPILE)dnl
dnl
process TS__PRC_NAME = TS__PRC_IN_PARTICLES => TS__PRC_OUT_PARTICLES
dnl
dnl Integration
dnl
divert(TS__DIV_INTEGRATE)dnl
dnl
printf "************************************************************************"
printf "* Integrating TS__PRC_DESCRIPTION"
printf "************************************************************************"
sqrts = TS__PRC_SQRTS
beams = TS__PRC_IN_PARTICLES
seed = 0
integrate (TS__PRC_NAME)
dnl
dnl Comparison
dnl
divert(TS__DIV_COMPARE)dnl
dnl
printf "************************************************************************"
printf "* Checking TS__PRC_DESCRIPTION"
printf "************************************************************************"
real error_sum = sqrt ((TS__PRC_ACCURACY) ** 2 + error(TS__PRC_NAME) ** 2)
printf "Expecting TS__PRC_INTEGRAL"
show(integral(TS__PRC_NAME))
real pull = abs (integral (TS__PRC_NAME) - TS__PRC_INTEGRAL) / error_sum
if (pull > 6) then
  printf "SEVERE:  pull > 6 in TS__PRC_DESCRIPTION"
elsif (pull > 5) then
  printf "SEVERE:  pull > 5 in TS__PRC_DESCRIPTION"
elsif (pull > 4) then
  printf "ERROR:   pull > 4 in TS__PRC_DESCRIPTION"
elsif (pull > 3) then
  printf "WARNING: pull > 3 in TS__PRC_DESCRIPTION"
elsif (pull > 2) then
  printf "NOTICE:  pull > 2 in TS__PRC_DESCRIPTION"
endif
tolerance = TS__PRC_TOLERANCE * error_sum
expect (integral (TS__PRC_NAME) == TS__PRC_INTEGRAL)
divert')dnl
dnl
dnl ====================================================================
dnl
dnl PARAMETERS (sindarin-statements-inserted-among-integration-steps)
dnl
dnl --------------------------------------------------------------------
define(`PARAMETERS',`dnl
dnl ====================================================================
divert(TS__DIV_INTEGRATE)dnl
$*
divert')dnl
dnl
dnl ====================================================================
dnl
dnl BEGIN_TESTSUITE (filename, [prefix])
dnl
dnl --------------------------------------------------------------------
define(`BEGIN_TESTSUITE',`dnl
dnl ====================================================================
dnl
ifelse(`$2',`',`',`define(`TS__PRC_PREFIX',`$2')')dnl
dnl
! Whizard test suite.  Do not edit.  Generated automatically from
! $1
! by the macros in
! $Id: testsuite.m4 5205 2014-01-30 14:17:39Z jr_reuter $
! ----------------------------------------------------------------------')
dnl ====================================================================
dnl
dnl END_TESTSUITE ()
dnl
dnl --------------------------------------------------------------------
define(`END_TESTSUITE',`dnl
dnl ====================================================================
! ----------------------------------------------------------------------
! Define the process
! ----------------------------------------------------------------------
undivert(TS__DIV_COMPILE)dnl
! ----------------------------------------------------------------------
! Compile the processes
! ----------------------------------------------------------------------
compile
! ----------------------------------------------------------------------
! Integrate the processes
! ----------------------------------------------------------------------
undivert(TS__DIV_INTEGRATE)dnl
! ----------------------------------------------------------------------
! Check the processes
! ----------------------------------------------------------------------
undivert(TS__DIV_COMPARE)dnl
! ----------------------------------------------------------------------
! Done
! ----------------------------------------------------------------------
printf "************************************************************************"
printf "* Done."
printf "************************************************************************"')
dnl
dnl ====================================================================
