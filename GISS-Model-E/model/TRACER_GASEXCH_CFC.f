      SUBROUTINE TRACERS_GASEXCH_ocean_CFC_PBL(tg1,ws,
     . alati,psurf,trmm,trconstflx,byrho,Kw_gas,alpha_gas,
     . beta_gas,trsf,trcnst,ilong,jlat)

      USE CONSTANT, only:    rhows,mair
      
      implicit none

      integer, intent(in) :: ilong,jlat
      real*8, intent(in)  :: trmm,tg1,ws,alati, psurf, trconstflx, byrho
      real*8, intent(out) :: Kw_gas, alpha_gas, beta_gas, trsf, trcnst
      real*8  :: Sc_gas
      real*8, external :: sc_cfc,sol_cfc
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !OCMIP implementation www.ipsl.jussieu.fr/OCMIP
      !F=Kw*Csat - Kw*Csurf=
      !  Kw*alpha*trs - Kw*trs
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !new treatment for units
      !F=Kw*Csat -Kw*Csurf
      ! =Kw*alpha*mol_weight_air/mol_weight_cfc11*surfp*Cair  !!!TERM_1*Cair
      ! -Kw*rho_water/mol_weight_cfc11*Csurf                  !!!TERM_2
      !
      ! where, Kw                in m/s
      !        alpha                mol/m^3/atm
      !        mol_weight_air       Kg_air
      !        mol_weight_cfc11     Kg_CFC-11
      !        surfp                atm
      !        Cair                 Kg_CFC-11/Kg_air
      !        rho_water            Kg_water/m^3
      !        Csurf                Kg_CFC-11/Kg_water
      !then F is in  (mol/m^3)(m/s)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !---------------------------------------------------------------
      !TRANSFER VELOCITY
      !---------------------------------------------------------------
      ! compute Schmidt number for gas
      ! use ground temperature in deg C including skin effects
      Sc_gas=sc_cfc(tg1,11)

      !wind speed ws: magn. of surf. wind modified by buoyancy flux (m/s)
      !compute transfer velocity Kw
      !only over ocean

      if (Sc_gas .le. 0.) then
        write(*,'(a,2i4,a,2f9.3)')
     .          'warning: Sc_gas negtv, at ',ilong,jlat,
     .          ', Sc_gas,temp_c=',Sc_gas,tg1
         Kw_gas=1.e-10
      else
         Kw_gas=
     &       1.d0/3.6e+5*0.337d0*ws*ws*(Sc_gas/660.d0)**(-0.5d0)
      endif


      !---------------------------------------------------------------
      !gas SOLUBILITY
      !---------------------------------------------------------------
      !alpha --solubility of CFC (11 or 12) in seawater
      !in mol/m^3/picoatm
       alpha_gas=sol_cfc(tg1,alati,11)
      !convert to mol/m^3/atm
       alpha_gas=alpha_gas*1.e+12

      !---------------------------------------------------------------
      !psurf is in mb. multiply with 10.197e-4 to get atm
      !include molecular weights for air and CFC-11
       beta_gas=alpha_gas*(psurf*10.197e-4)*mair*1.e-3
     &                   /(trmm*1.e-3)
!!!    beta_gas = beta_gas * trmm*1.e-3/rhows

      !trsf is really sfac = Kw_gas * beta_gas
      !units are such that flux comes out to (m/s)(kg/kg)
       trsf = Kw_gas * beta_gas

       trcnst = Kw_gas * trconstflx*byrho ! convert to (conc * m/s)

      RETURN
      END SUBROUTINE TRACERS_GASEXCH_ocean_CFC_PBL

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------

c used with TRACERS_GASEXCH_ocean_CFC to compute transfer velocity for CFCs
c
c ---------------------------------------------------------------------
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
      a1 ( 11) = -229.9261d0
      a2 ( 11) =  319.6552d0
      a3 ( 11) =  119.4471d0
      a4 ( 11) =   -1.39165d0
      b1 ( 11) =   -0.142382d0
      b2 ( 11) =    0.091459d0
      b3 ( 11) =   -0.0157274d0
c
c     for CFC/12
c     ----------
      a1 ( 12) = -218.0971d0
      a2 ( 12) =  298.9702d0
      a3 ( 12) =  113.8049d0
      a4 ( 12) =   -1.39165d0
      b1 ( 12) =   -0.143566d0
      b2 ( 12) =    0.091015d0
      b3 ( 12) =   -0.0153924d0
c
      ta       = ( pt + tf)* 0.01d0
      d    = ( b3 ( kn)* ta + b2 ( kn))* ta + b1 ( kn)

c
c
      sol_cfc
     $    = exp ( a1 ( kn)
     $    +       a2 ( kn)/ ta
     $    +       a3 ( kn)* dlog ( ta )
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

      REAL*8 FUNCTION alpha_gas2_cfc(pt,ps)
!@sum helper function for SURFACE calculation
      real*8 :: pt,ps,sol_cfc

      alpha_gas2_cfc=1e12*sol_cfc(pt,ps,11) !convert to mol/m^3/atm

      return
      end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c used with TRACERS_GASEXCH_ocean_CFC to compute transfer velocity for CFCs
c
c  ---------------------------------------------------------------------
      REAL*8 FUNCTION sc_cfc(t,kn)
c---------------------------------------------------
c     CFC 11 and 12 Schmidt number 
c     as a function of temperature. 
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
c   coefficients with t in degre Celcius
c   ------------------------------------
      a1(11) = 3501.8d0
      a2(11) = -210.31d0
      a3(11) =    6.1851d0
      a4(11) =   -0.07513d0
c
      a1(12) = 3845.4d0
      a2(12) = -228.95d0
      a3(12) =    6.1908d0
      a4(12) =   -0.067430d0
c

      sc_cfc = a1(kn) + a2(kn) * t + a3(kn) *t*t  
     &         + a4(kn) *t*t*t
  
      RETURN 
      END 
