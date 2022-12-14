#include "rundeck_opts.h"
      subroutine OCN_cfcflx
c
      USE MODEL_COM, only: dtsrc
      USE obio_forc, only: wind,atmCFC
      use TimeConstants_mod, only: SECONDS_PER_HOUR, DAYS_PER_YEAR,
     &                             HOURS_PER_DAY
      
     .                    ,ao_cfcflux

      USE MODEL_COM, only : nstep=>itime
      USE OCEANRES, only : kdm=>lmo

      implicit none


      integer :: i,j,k

      integer :: nt,kmax
      real  :: sccfc,sccfcarg,wssq,rkwcfc
      real  :: Ts,tk,ff,xcfc,deltcfc,flxmolm3,flxmolm3h

!this is for gas exchange done in the ocean grid: 
      k = 1
      Ts = temp1d(k)
      Ss = saln1d(k)
      sccfc = sc_cfc(Ts,11)   !for CFC-11   %from OCMIP-2

      wssq = wind*wind
      if (sccfc.lt.0.) then
        sccfcarg=1.d-10
        rkwcfc=1.d-10
      else
        sccfcarg = (sccfc/660.D0)**(-0.5)     !Schmidt number
        awan=0.337/(3.6E+5)                   !piston vel coeff., from
                                              !Wanninkof 1992, but adjusted
                                              !by OCMIP, and converted from
                                              !cm/hr to m/s
        rkwcfc = awan*wssq*sccfcarg           !transfer coeff (units of m/s)
      endif
      ff = sol_cfc(Ts,Ss,11)      !for CFC-11, in units mol/m3/pptv =mol/m3/uatm
                                        !parts per trillion by volume (1 pptv  =  mole fraction 1 x 10-12)

      atmCFC = 263.                     ! for pCFC-11, year 1998, units pptv (uatm)
      xcfc = atmCFC*1013.D0/stdslp
      !make sure pCFC is also in uatm
      deltcfc = (xcfc-pCFC_ij)*ff           !units in mol/m3
      flxmolm3 = 1000.d0 * (rkwcfc*deltcfc/dp1d(k))   !units of mili-mol/m3/s
      term = flxmolm3*pnoice(k)    !units of uM/s (=mili-mol/m^3/s)

      !flux sign is (atmos-ocean)>0, i.e. positive flux is INTO the ocean
      ao_cfcflux= rkwcfc*(xcfc-pCFC_ij)*ff*pnoice(k)            ! air-sea cfc flux
     .            *SECONDS_PER_HOUR                             ! mol/m2/hr
     .            *137.37*HOURS_PER_DAY*DAYS_PER_YEAR           ! grC/m2/yr


      return
      end

c ---------------------------------------------------------------------------
c_ $Source: /home/orr/ocmip/web/OCMIP/phase2/simulations/CFC/boundcond/RCS/sc_cfc.f,v $
c_ $Revision: 1.1 $
c_ $Date: 1998/07/07 15:22:00 $   ;  $State: Exp $
c_ $Author: orr $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: sc_cfc.f,v $
c_ Revision 1.1  1998/07/07 15:22:00  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_
      REAL*8 FUNCTION sc_cfc(t,kn)
c---------------------------------------------------
c     CFC 11 and 12 schmidt number 
c     as a fonction of temperature. 
c
c     ref: Zheng et al (1998), JGR, vol 103,No C1 
c
c     t: temperature (degree Celcius)
c     kn: = 11 for CFC-11,  12 for CFC-12
c
c     J-C Dutay - LSCE
c---------------------------------------------------
      IMPLICIT NONE 
      INTEGER kn
      REAL*8  a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      REAL*8 t
c
c   coefficients with t in degree Celcius
c   ------------------------------------
      a1(11) = 3501.8
      a2(11) = -210.31
      a3(11) =    6.1851
      a4(11) =   -0.07513
c
      a1(12) = 3845.4
      a2(12) = -228.95
      a3(12) =    6.1908
      a4(12) =   -0.067430
c

      sc_cfc = a1(kn) + a2(kn) * t + a3(kn) *t*t  
     &         + a4(kn) *t*t*t
  
      RETURN 
      END 
c ---------------------------------------------------------------------------
c_ $Source: /www/html/ipsl/OCMIP/phase2/simulations/CFC/boundcond/RCS/sol_cfc.f,v $
c_ $Revision: 1.2 $
c_ $Date: 1998/07/17 07:37:02 $   ;  $State: Exp $
c_ $Author: jomce $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: sol_cfc.f,v $
c_ Revision 1.2  1998/07/17 07:37:02  jomce
c_ Fixed slight bug in units conversion: converted 1.0*e-12 to 1.0e-12
c_ following warning from Matthew Hecht at NCAR.
c_
c_ Revision 1.1  1998/07/07 15:22:00  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_
      REAL*8 FUNCTION sol_cfc(pt,ps,kn)
c-------------------------------------------------------------------
c
c     CFC 11 and 12 Solubilities in seawater
c     ref: Warner & Weiss (1985) , Deep Sea Research, vol32
c
c     pt:       temperature (degre Celcius)
c     ps:       salinity    (o/oo)
c     kn:       11 = CFC-11, 12 = CFC-12
c     sol_cfc:  in mol/m3/pptv
c               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm

c
c     J-C Dutay - LSCE
c-------------------------------------------------------------------
      USE CONSTANT, only : tf

      REAL*8    pt, ps,ta,d
      REAL*8 a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      REAL*8 b1 ( 11: 12), b2 ( 11: 12), b3 ( 11: 12)

      INTEGER kn

cc
cc coefficient for solubility in  mol/l/atm
cc ----------------------------------------
c
c     for CFC 11
c     ----------
      a1 ( 11) = -229.9261
      a2 ( 11) =  319.6552
      a3 ( 11) =  119.4471
      a4 ( 11) =   -1.39165
      b1 ( 11) =   -0.142382
      b2 ( 11) =    0.091459
      b3 ( 11) =   -0.0157274
c    
c     for CFC/12
c     ----------
      a1 ( 12) = -218.0971
      a2 ( 12) =  298.9702
      a3 ( 12) =  113.8049
      a4 ( 12) =   -1.39165
      b1 ( 12) =   -0.143566
      b2 ( 12) =    0.091015
      b3 ( 12) =   -0.0153924
c

      ta       = ( pt + tf)* 0.01
      d    = ( b3 ( kn)* ta + b2 ( kn))* ta + b1 ( kn)
c
c
      sol_cfc 
     $    = exp ( a1 ( kn)
     $    +       a2 ( kn)/ ta 
     $    +       a3 ( kn)* alog ( ta )
     $    +       a4 ( kn)* ta * ta  + ps* d )
c
c     conversion from mol/(l * atm) to mol/(m^3 * atm) 
c     ------------------------------------------------
      sol_cfc = 1000. * sol_cfc
c
c     conversion from mol/(m^3 * atm) to mol/(m3 * pptv) 
c     --------------------------------------------------
      sol_cfc = 1.0e-12 * sol_cfc

      END 
      
