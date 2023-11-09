#include "rundeck_opts.h"

module tracers_dCO

implicit none

!@param dacetone_fact factor to multiply with acetone (scaled) conc
!@param dalke_IC_fact factor to multiply with Alkenes initial conc
!@param dPAR_IC_fact factor to multiply with Paraffin initial conc
!@param dPAN_IC_fact factor to multiply with PAN initial conc
!@param dMeOOH_IC_fact factor to multiply with MeOOH initial conc
!@param dHCHO_IC_fact factor to multiply with HCHO initial conc
!@param dCO_IC_fact factor to multiply with CO_IC file
!@param d17O2_to_O2 Ratio of d17O2 over total O2 in the atmosphere.
!@param d18O2_to_O2 Ratio of d18O2 over total O2 in the atmosphere.

#ifdef TRACERS_dCO_bin_reprod
real*8, parameter :: dacetone_fact=1.d0
real*8, parameter :: dalke_IC_fact=1.d0
real*8, parameter :: dPAR_IC_fact=1.d0
real*8, parameter :: dPAN_IC_fact=1.d0
real*8, parameter :: dMeOOH_IC_fact=1.d0
real*8, parameter :: dHCHO_IC_fact=1.d0
real*8, parameter :: dCO_IC_fact=1.d0
real*8, parameter :: d17O2_to_O2=1.d0
real*8, parameter :: d18O2_to_O2=1.d0
#else
you should not even compile, not implemented yet.
#endif  /* TRACERS_dCO_bin_reprod */

end module tracers_dCO
