#include "rundeck_opts.h"

      SUBROUTINE TRACERS_GASEXCH_ocean_CO2_PBL(tg1,ws,
     . alati,psurf,trmm,trconstflx,byrho,Kw_gas,alpha_gas,
     . beta_gas,trsf,trcnst,ilong,jlat)

      USE CONSTANT, only:    rhows
      use OldTracer_mod, only: vol2mass

      USE TRACER_COM, only : gasex_index, n_co2n
      
      implicit none

      real*8, parameter :: awan=0.337d0/(3.6d5) !piston vel coeff., from
                                                !Wanninkof 1992, but adjusted
                                                !by OCMIP, and converted from
                                                !cm/hr to m/s
      integer, intent(in) :: ilong,jlat
      real*8, intent(in)  :: trmm,tg1,ws,alati, psurf, trconstflx, byrho
      real*8, intent(out) :: Kw_gas, alpha_gas, beta_gas, trsf, trcnst
      real*8  :: Sc_gas
      real*8, external :: sc_co2,sol_co2
      integer :: idx

!@var  psurf surface pressure
!@var  alati SSS at i,j
!@var  Kw_gas  gas exchange transfer velocity at i,j only over ocean
!@var  alpha_gas  solubility of gas
!@var  beta_gas  conversion term  that includes solubility
!@var  trconstflx  constant component of surface tracer flux
!      byrho   1/surface air density

!routine for PBL calculations of CO2 gas exchange

      idx=gasex_index%getindex(n_co2n)
      !---------------------------------------------------------------
      !TRANSFER VELOCITY
      !---------------------------------------------------------------
      !Schmidt number for gas
       Sc_gas=sc_co2(tg1)

             !wind speed ws: magn. of surf. wind modified by buoyancy flux (m/s)
      !compute transfer velocity Kw only over ocean
      if (Sc_gas .le. 0.) then
!        write(*,'(a,2i4,a,2f9.3)')
!     .          'warning: Sc_gas negtv, at ',ilong,jlat,
!     .          ', Sc_gas,temp_c=',Sc_gas,tg1
         Kw_gas=1.e-10
      else
         Kw_gas=(Sc_gas/660.d0)**(-0.5d0) * ws * ws * awan !units of m/s
      endif

      !---------------------------------------------------------------
      !gas SOLUBILITY
      !---------------------------------------------------------------
      !alpha --solubility of CO2 in seawater
      alpha_gas=sol_co2(tg1,alati)    !mol/m^3/picoatm
      alpha_gas = alpha_gas * 1.d-6   !mol,CO2/m3/uatm
      alpha_gas = alpha_gas * 1024.5  !Watson Gregg
      !---------------------------------------------------------------
      !psurf is in mb. 
      beta_gas = alpha_gas * psurf/1013.25      !stdslp and psurf in mb, no need to change units

      !trsf is really sfac = Kw_gas * beta_gas
      !the term 1.0d6 / vol2mass(idx) is needed to convert uatm -> kg,co2/kg,air 
      !in the denominator of alpha
      trsf = Kw_gas * beta_gas * 1.0d6 / vol2mass(idx)

      !trconstflx comes in from SURFACE.f and has units kg,co2/kg,air/m2
      !therefore trcnst needs to be multiplied by byrho before it is sent to  PBL.f
      trcnst = Kw_gas * alpha_gas * trconstflx * byrho     
     .                * 1.0d6 / vol2mass(idx)    

#ifndef OBIO_QUIET_MODE
        if (ilong.eq.1. .and. jlat.eq.45) then
       write(*,'(a,2i7,11e12.4)')'PBL, TRACER_GASEXCH_CO2 ws:',
!       write(*,'(a,2i7,11e12.4)')'44444444444444444444444',  
     .   ilong,jlat,tg1,(Sc_gas/660.d0)**(-0.5d0),ws*ws,
     .   Kw_gas,alpha_gas,beta_gas,trsf,trcnst,trconstflx,byrho,rhows
        endif
#endif /* not OBIO_QUIET_MODE */

      RETURN
      END SUBROUTINE TRACERS_GASEXCH_ocean_CO2_PBL

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c used with TRACERS_GASEXCH_ocean_CO2 to compute transfer velocity for CO2
c
c taken from Watsons' code
c ---------------------------------------------------------------------
      REAL*8 FUNCTION sol_co2(pt,ps)
c-------------------------------------------------------------------
c
c     CO2 Solubility in seawater
c
c     pt:       temperature (degrees Celcius)
c     ps:       salinity    (o/oo)
c     sol_co2:  in mol/m3/pptv
c               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm
c-------------------------------------------------------------------
      USE CONSTANT, only : tf
      REAL*8    pt,ps,ta,tk100,tk1002


      ta = (pt + tf)
      tk100 = ta*0.01d0
      tk1002 = tk100*tk100
!     sol_co2 = exp(-162.8301d0 + 218.2968d0/tk100  +
!    .          90.9241d0*log(tk100) - 1.47696d0*tk1002 +
!    .          ps * (.025695d0 - .025225d0*tk100 +
!    .          0.0049867d0*tk1002))
      !new OCMIP2016 values 
      sol_co2 = exp(-160.7333d0 + 215.4152d0/tk100  +
     .          89.8920d0*log(tk100) - 1.47759d0*tk1002 +
     .          ps * (0.029941d0 - 0.027455d0*tk100 +
     .          0.0053407d0*tk1002))

      END

      REAL*8 FUNCTION alpha_gas2_co2(pt,ps)
c-------------------------------------------------------------------
c
c     CO2 Solubility in seawater
c
c     pt:       temperature (degrees Celcius)
c     ps:       salinity    (o/oo)
c     sol_co2:  in mol/m3/pptv
c               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm
c     Scaled by 1024.5 (density of sea water) and 10^6 to get
c     alpha_gas in  mol/m3/uatm
c-------------------------------------------------------------------

      REAL*8    pt,ps,sol_CO2

      alpha_gas2_co2=sol_CO2(pt,ps) * 1.d-6 * 1024.5  !mol,CO2/m3/uatm
      
      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      REAL*8 FUNCTION sc_co2(t)
c---------------------------------------------------
c     CO2 Schmidt number 
c     as a function of temperature. 
c
c     t: temperature (degree Celcius)
c---------------------------------------------------
      IMPLICIT NONE 
      REAL*8  a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      REAL*8 t
!
!     sc_co2 = 2073.1 - 125.62*t + 3.6276*t**2 - 0.043219*t**3
      !new OCMIP2016 values from Wanninkhof (2014)
      sc_co2 = 2116.8d0 - 136.25d0*t + 4.7353d0*t**2 - 0.092307d0*t**3
     &          + 0.0007555d0*t**4
!
      RETURN 
      END 

