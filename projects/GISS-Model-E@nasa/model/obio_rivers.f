#include "rundeck_opts.h"

#ifdef OBIO_RUNOFF
      subroutine obio_rivers(vrbos)

      use MODEL_COM, only: dtsrc,nstep=>itime
      use obio_incom, only: estFe
      use obio_com, only:  vrbos,rpocconc,rnitrconc,rsiliconc
     .                    ,rironconc,rdocconc,rdicconc
!                         ,rnitrfmlo
     .                    ,rhs,P_tend,D_tend,C_tend,dp1d
#ifdef TRACERS_Alkalinity
     .                    ,ralkconc,A_tend
#endif
      use obio_forc, only: river_runoff

      implicit none
    
      real :: term
      logical :: vrbos

         term = rnitrconc
     .    * river_runoff/dtsrc        ! kg,N/kg,w => kg,N/m2,w/s
     .    / dp1d(1)                   ! kg,N/m2,w/s => kg,N/m3,w/s
     .    * 1.d6/14.                  ! kg,N/m3,w/s => mmol,N/m3,w/s
         rhs(1,1,17) = term
         P_tend(1,1) = P_tend(1,1) + term

        term = rsiliconc
     .   * river_runoff/dtsrc        ! kg,S/kg,w => kg,S/m2,w/s
     .   / dp1d(1)                   ! kg,S/m2,w/s => kg,S/m3,w/s
     .   * 1.d6/28.055               ! kg,S/m3,w/s => mmol,S/m3,w/s
        rhs(1,3,17) = term
        P_tend(1,3) = P_tend(1,3) + term

        term = rironconc
     .   * river_runoff/dtsrc        ! kg,Fe/kg,w => kg,Fe/m2,w/s
     .   / dp1d(1)                   ! kg,Fe/m2,w/s => kg,Fe/m3,w/s
     .   * 1.d9/55.845               ! kg,Fe/m3,w/s => umol,Fe/m3,w/s
     .   * estFe                     ! estuarine retention rate
        rhs(1,4,17) = term
        P_tend(1,4) = P_tend(1,4) + term

        term = rpocconc
     .   * river_runoff/dtsrc        ! kg,C/kg,w => kg,C/m2,w/s
     .   / dp1d(1)                   ! kg,C/m2,w/s => kg,C/m3,w/s
     .   * 1.d6                      ! kg,C/m3,w/s => mg,C/m3,w/s
        rhs(1,10,17) = term
        D_tend(1,1) = D_tend(1,1) + term

        term = rdocconc
     .    * river_runoff/dtsrc         ! kg,C/kg,w => kg,C/m2,w/s
     .    / dp1d(1)                    ! kg,C/m2,w/s => kg,C/m3,w/s
     .    * 1.d6/12.                   ! kg,C/m3,w/s => mmol,C/m3,w/s
        rhs(1,13,17) = term
        C_tend(1,1) = C_tend(1,1) + term

        term = rdicconc
     .    * river_runoff/dtsrc         ! kg,C/kg,w => kg,C/m2,w/s
     .    / dp1d(1)                    ! kg,C/m2,w/s => kg,C/m3,w/s
     .    * 1.d6/12.                   ! kg,C/m3,w/s => mmol,C/m3,w/s
        rhs(1,14,17) = term
        C_tend(1,2) = C_tend(1,2) + term

#ifdef TRACERS_Alkalinity
      term = ralkconc
     . * river_runoff/dtsrc          ! mol/kg,w => mol/m2,w/s
     . / dp1d(1)                     ! mol/m2,w/s => mol/m3,w/s
     . * 1.d3                        ! mol/m3,w/s => mmol/m3,w/s
      rhs(1,15,17) = term
      A_tend(1) = A_tend(1) + term
#endif 

      end subroutine obio_rivers
#endif
