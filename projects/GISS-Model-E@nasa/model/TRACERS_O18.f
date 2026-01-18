#include "rundeck_opts.h"

!@sum  TRACERS_O18 water isotope specific routines and functions
!@auth Gavin Schmidt

      FUNCTION FRACVL(TEMP,itr)
!@sum FRACVL Calculate vapor-->liquid equilibrium fractionation factor
!@auth Gavin Schmidt
      USE CONSTANT, only : tf
      use oldtracer_mod, only : iso_index
      IMPLICIT NONE
!@var TEMP  temperature (deg C)
      REAL*8, INTENT(IN) :: TEMP
      INTEGER, INTENT(IN) :: itr  ! actual tracer number
      INTEGER, PARAMETER :: NTSPM=5
      REAL*8, PARAMETER ::
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium 5: o17
c**** exponential fit 
     * A(NTSPM)=(/  0d0,  1137d0   ,  24844d0  , 46480d0,   601.47d0/),
     * B(NTSPM)=(/  0d0, -0.4156d0 , -76.248d0 ,-103.87d0, -0.2199d0/),
     * C(NTSPM)=(/  0d0, -2.0667d-3,  52.612d-3, 0d0,      -1.0933d-3/)
c**** quadratic fit to Majoube (1971)
c     *     A(NTSPM) = (/  0d0, -3.75d-7, -6.375d-6, -6.875d-6  /),
c     *     B(NTSPM) = (/  0d0,  1.025d-4, 1.2475d-3, 1.7d-3    /),
c     *     C(NTSPM) = (/  1d0,  0.9884d0, 0.9001d0 , 0.86975d0 /)
      REAL*8 FRACVL,TK
      INTEGER N
C****
      N=iso_index(itr)
C**** Quadratic fit
c      FRACVL=C(N) + TEMP*(B(N) + TEMP*A(N))
C**** Exponential
      TK=TEMP+TF
      FRACVL=EXP(-A(N)/TK**2 - B(N)/TK - C(N))
C****
      RETURN
      END FUNCTION FRACVL

      FUNCTION FRACVS(TEMP,itr)
!@sum FRACVS Calculate vapour --> solid (ice) equil. fractionation fact.
!@auth Gavin Schmidt
      USE CONSTANT, only : tf
      use oldtracer_mod, only : iso_index
      IMPLICIT NONE
!@var TEMP  temperature (deg C)
      REAL*8, INTENT(IN) :: TEMP
      INTEGER, INTENT(IN) :: itr ! actual tracer number
      INTEGER, PARAMETER :: NTSPM=5
      REAL*8, PARAMETER ::
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium 5: o17
C**** exponential
     * A(NTSPM) = (/ 0d0,  0d0       , 16288d0 , 46480d0 , 0d0/),
     * B(NTSPM) = (/ 0d0,  11.839d0  , 0d0     ,-103.87d0, 6.2628d0/),
     * C(NTSPM) = (/ 0d0, -0.028244d0,-0.0934d0, 0d0,     -.01494d0/)
C**** linear fit
c     *     A(NTSPM) = (/  0d0,  1.36d-4,   1.46d-3, 0d0/),
c     *     B(NTSPM) = (/  1d0,  0.9850d0, 0.8834d0, 0.7845d0/)

      REAL*8 FRACVS,TK
      INTEGER N

      N=iso_index(itr)
C**** linear fit
c      FRACVS=B(N) + A(N)*TEMP
C**** Exponential
      TK=TEMP+TF
      FRACVS=EXP(-A(N)/TK**2 - B(N)/TK -C(N))
C****
      RETURN
      END FUNCTION FRACVS

      FUNCTION FRACLS(itr)
!@sum FRACLS Calculate liquid --> solid equilibrium fractionation factor
!@auth Gavin Schmidt
      use oldtracer_mod, only : iso_index
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: itr ! actual tracer number
      INTEGER, PARAMETER :: NTSPM=5
      REAL*8, PARAMETER ::
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium 5: o17
     *     A(NTSPM) = (/ 1d0, 1.0035d0, 1.0208d0,  1.04d0, 1.00185d0/)

      REAL*8 FRACLS
C****
      FRACLS=A(iso_index(itr))
C****
      RETURN
      END FUNCTION FRACLS

      FUNCTION FRACLK(WS,itr)
!@sum FRACLK calculates the liquid/vapor kinetic fractionation factor
!@+          from either (Merlivat and Jouzel,1979) or as a constant
!@auth Gavin Schmidt
      use oldtracer_mod, only : iso_index
      IMPLICIT NONE
!@var WS surface wind speed (m/s)
      REAL*8, INTENT(IN) :: WS
      INTEGER, INTENT(IN) :: itr  ! actual tracer index
      INTEGER, PARAMETER :: NTSPM=5
      REAL*8, PARAMETER ::
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium 5: o17
     * A(NTSPM)=(/ 1d0, 0.994d0,  0.99472d0,  0.98944d0,  0.996913d0 /),
     * B(NTSPM)=(/ 0d0, 0.285d-3, 0.2508d-3,  0.5016d-3,  0.14663d-3 /),
     * C(NTSPM)=(/ 0d0, 0.82d-3,  0.7216d-3,  0.14432d-2, 0.42189d-3 /)
C**** alternative from Cappa et al (2003)
c     * A(NTSPM) = (/ 1d0, 0.99295d0,  0.99636d0,  0.99295d0, 0.99637d0/),
c     * B(NTSPM) = (/ 0d0, 0.485d-3, 0.188d-3,  0.485d-3, .2495d-3/),
c     * C(NTSPM) = (/ 0d0, 0.727d-3, 0.275d-3,  0.727d-3, .3740d-3/)
      REAL*8 FRACLK
      INTEGER N
C****
C**** Calculate kinetic fractionation factor:
C****
      N=iso_index(itr)
      FRACLK=A(N)
      IF(WS.GE.7.) FRACLK=1d0-(B(N)*WS+C(N))
C****
      RETURN
      END FUNCTION FRACLK

#ifdef TRACERS_SPECIAL_O18
      SUBROUTINE ISOEQUIL(N,TEMP,LIQU,QMV,QML,TRMV,TRML,FEQ)
!@sum ISOEQUIL equilibrates isotopes in vapor & liquid/solid reservoirs
!@auth Gavin Schmidt/Georg Hoffmann
      USE CONSTANT, only : tf
      USE TRACER_COM, only: tr_wd_TYPE,nWATER
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: LIQU !@var LIQU true if QML is liquid
      REAL*8, INTENT(IN) :: TEMP  !@var TEMP temperature (K)
!@var QMV,QML vapour and liquid water masses (kg or kg/m^2)
      REAL*8, INTENT(IN) :: QMV,QML
!@var TRMV,TRML vapour and liquid tracer mass (kg or kg/m^2)
      REAL*8, INTENT(INOUT) :: TRMV,TRML
!@var FEQ fraction of condensate equilibrated (1. = all, 0. = none)
      REAL*8, INTENT(IN) :: FEQ
      REAL*8 TDEGC,ZALPH,ZDELEQU,ZXFAC,FRACVL,FRACVS

      SELECT CASE(tr_wd_TYPE(N))
      CASE(nWATER)
        TDEGC = TEMP - TF
        IF (QMV.GT.0.) THEN
          IF (LIQU) THEN
            ZALPH = 1./FRACVL(TDEGC,N)
          ELSE
            ZALPH = 1./FRACVS(TDEGC,N)
          END IF
          ZXFAC = FEQ*ZALPH*QML/QMV
          ZDELEQU = (FEQ*TRML - ZXFAC*TRMV)/(1.+ZXFAC)
          TRML = TRML - ZDELEQU
          TRMV = TRMV + ZDELEQU
        END IF
      END SELECT

      RETURN
C****
      END SUBROUTINE ISOEQUIL

      FUNCTION KIN_COND_ICE(ALPH,SUPSAT,itr)
!@sum calculate kinetic fractionation when condensing to ice in 
!@+   super-saturated conditions
!@auth Gavin Schmidt/Georg Hoffmann
      USE CONSTANT, only : tf
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: itr   ! actual tracer index
!@var SUPSAT super_saturation factor from cloud scheme
      REAL*8, INTENT(IN) :: SUPSAT
!@var ALPH equilibrium fractionation
      REAL*8, INTENT(IN) :: ALPH
      real*8 kin_cond_ice, get_diff_rel
C****
C**** Calculate kinetic condensation when condensing to ice
C****
      KIN_COND_ICE=alph*SUPSAT/(1.+(SUPSAT-1.)*alph*get_diff_rel(itr))

      return
      end function kin_cond_ice

      FUNCTION KIN_EVAP_PREC(ALPH,HEFF,itr)
!@sum calculate kinetic fractionation when evaporating into 
!@+   undersaturated environment
!@auth Gavin Schmidt/Georg Hoffmann, 
!@+    Stewart/Rayleigh limit added by Jesse Nusbaumer
      USE CONSTANT, only : tf
      use oldtracer_mod, only : iso_index
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: itr ! actual tracer index
!@var HEFF effective relative humidity from cloud scheme
      REAL*8, INTENT(IN) :: HEFF
!@var alph equilibrium fractionation 
      REAL*8, INTENT(IN) :: ALPH
      INTEGER, PARAMETER :: NTSPM=5
!@var ZDIFRELGAM = ZDIFREL^0.58 
      real*8 :: ZDIFRELGAM(NTSPM) =
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium 5: o17
     *     (/ 1d0, 1.0164d0, 1.0145d0, 1.0191d0, 1.008479d0 /)
C**** alternative from Cappa et al (2003)
c     *     (/ 1d0, 1.0184d0, 1.0095d0, 1.0184d0, 1.009478d0 /)
!@var raylimit Stewart/Rayleigh fractionation factor at HEFF=0
      real*8 raylimit
!
!
      real*8 kin_evap_prec
C****
C**** Calculate kinetic condensation when evaporating below clouds
C****
      KIN_EVAP_PREC=alph*HEFF/(1.+(HEFF-1.)*alph
     *     *ZDIFRELGAM(iso_index(ITR)))
      
!     ------------------------------------------------------------
!     The equation above produes questioniable and/or non-physical
!     results at low relative humidities. This may be due to the 
!     the fact that the equation was originally derived for
!     vapor depositing onto ice in super-saturated conditions,
!     not rain evaporatiing into low relative humidity conditions,
!     and thus some of the assumptions used to derive the equation
!     may be invalid.  To avoid these "bad" values, it is assumed
!     that the largest distillation effect possible is with a pure
!     Rayleigh distillation as derived by Stewart, 1975 for the limit
!     when the relative humidity is exactly zero.
!     ------------------------------------------------------------

      !Calculate the Rayleigh limit (e.g. HEFF=0) as derived
      !by Stewart, 1975
       raylimit = alph/ZDIFRELGAM(iso_index(itr))

      !If original fractionation value is less than
      !(i.e. more fractionating than) Rayleigh limit,
      !set back to Rayleigh limit:
       if(kin_evap_prec .lt. raylimit) kin_evap_prec = raylimit

      !If original fractionation value is greater than or
      !equal to one (which implies the droplet became more
      !depleted during distillation, which is unphysical),
      !replace it with the Rayleigh limit value:
       if(kin_evap_prec .ge. 1.d0) kin_evap_prec = raylimit  

!     ------------------------------------------------------------  

      return
      end function kin_evap_prec

      FUNCTION delta(W,T,n)
!@sum calculate per mil values
      use TRACER_COM, only : trw0
      implicit none
      real*8, intent(in) :: W,T
      integer, intent(in) :: n
      real*8 delta

      delta=0.
      if (trw0(n).gt.0 .and. W.gt.0) delta=(T/(W*trw0(n))-1.)*1d3

      return
      end function delta

      subroutine get_frac(itype,ws,tg1,trcnst,trsf,evap,tgr,q1,itr
     $     ,fac_cq_tr,trc1,trs1,fk)
!@sum get_frac calculate fractionation factors during evaporation
!@auth Gavin Schmidt
      use constant, only : tf
      implicit none
!@var fac_cq_tr kinetic frac. explicit Schmidt no. dependence
      real*8, intent(in) :: ws,tg1,tgr,trcnst,trsf,fac_cq_tr,evap,q1
      integer, intent(in) :: itype
      integer, intent(in) :: itr ! actual tracer number
!@var trc1 factor multiplying trcnst in PBL
!@var trs1 factor multiplying trsf in PBL
      real*8, intent(out) :: trc1,trs1
      real*8 frac,fk,fracvl,fracvs,fraclk

C**** Isotope tracers have different fractionations dependent on
C**** type and direction of flux
      select case (itype)
      case (1)                  ! ocean: kinetic fractionation
#ifdef O18_KINETIC_FRAC
        fk = fac_cq_tr
#else
        fk = fraclk(ws,itr)
#endif
        trc1 = trcnst * fk * fracvl(tg1,itr)
        trs1 = trsf * fk
      case (2:3)              ! other types
C**** tracers are now passive, so use 'upstream' concentration
        if (evap.lt.0) then ! dew
          trc1 = 0.
          if (tg1.gt.0) then
            frac=fracvl(tg1,itr)
          else
            frac=fracvs(tg1,itr)
          end if
          trs1=-evap/(q1*frac)
        else
          trc1 = evap*tgr
          trs1 = 0.
        end if
        fk=1.0  ! diagnostic only
      case (4)
C**** land surface tracers have correct values already
        trc1=trcnst
        trs1=trsf
      end select
      return
      end subroutine get_frac

      real*8 function get_diff_rel(itr)
      use oldtracer_mod, only : iso_index
      implicit none
      integer, intent(in) :: itr  ! tracer index

      integer, parameter :: ntspm=5
!@var zdifrel = inverse ratio of diffusion coeffs w.r.t normal water
      real*8 :: zdifrel(ntspm) =
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium 5: o17
! MJ78
     *     (/ 1d0 ,1.0285d0, 1.0251d0, 1.0331d0, 1.014663d0/) 
! Cappa et al 2003
c     *     (/ 1d0 ,1.0319d0, 1.0164d0, 1.0319d0, 1.016399d0/) 

      get_diff_rel = ZDIFREL(iso_index(itr))

      return
      end function
#endif
      
      function calc_permil(Rsam,Msam,nTr)
c**** Calculates concentration in per mille units
      
      USE TRACER_COM, only : trw0 ! amount in standard

      real*8 Rsam ! amount in sample
     *     ,Msam  ! mass of sample (amount of 'water' tracer)
      integer nTr ! sample tracer number

      calc_permil = 1.d3*(Rsam/(Msam*trw0(nTr))-1.)

      return
      end function
     
      function water_iso_conc_in_aquifer(n) Result(c)
!@sum returns concentration of tracer n in aquifer (kg/kg)
!@+   This is a slow function. Do not use it in loops!
      USE TRACER_COM, only : trw0
      use OldTracer_mod, only: trname
      USE Dictionary_mod, only : get_param
      implicit none
      integer, intent(in) :: n !@var tracer number
      real*8 :: c !@var concentration of this tracer in water (kg/kg)
      real*8 :: in_permil

      select case(trname(n))
      case('H2O18')
        call get_param("H2O18_in_aquifer", in_permil, default=-8.d0)
      case('HDO')
        call get_param("HDO_in_aquifer", in_permil, default=-54.d0)
      case('H2O17')
        !Default assumes 17O-excess of zero:
        call get_param("H2O17_in_aquifer", in_permil, 
     &                 default=-4.240036d0)
      case default
        c = 0.d0
        return
      end select

      c = (1.d0 + in_permil*1.d-3) * trw0(n)

      end function water_iso_conc_in_aquifer

