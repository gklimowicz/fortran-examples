C****  
C**** OCNKPP.f    KPP OCeaN vertical mixing scheme    2009/08/25
C****
#include "rundeck_opts.h"

      MODULE KPP_COM
!@sum  KPP_COM holds variables related to the KPP mixing scheme
!@auth Gavin Schmidt
      USE OCEAN, only : im,jm,lmo,kpl
      USE SW2OCEAN, only : lsrpd
      IMPLICIT NONE
      SAVE

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::    G0M1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  MO1,GXM1,SXM1, UO1,UOD1
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: S0M1,GYM1,SYM1, VO1,VOD1
#ifdef TRACERS_OCEAN
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRMO1,TXMO1,TYMO1
#endif

!@dbparam use_tdiss whether to apply tidally induced vertical mixing
!@+       from input dataset
      integer :: use_tdiss=0
!@dbparam tdiss_eff locality+conversion efficiency of dissip -> kv
      real*8 :: tdiss_eff=.7d0/3d0

!@var tdiss prescribed tidal dissipation (W/m2)
!@var tdiss_n N corresponding to prescribed tidal dissipation (1/s)
      real*8, allocatable, dimension(:,:) :: tdiss,tdiss_n

      ! todo: make fkph a runtime parameter along with difsiw in KPPE
      real*8, parameter ::
     &    fkph   = 0.1d0    ! background vertical diffusion coefficient (cm**2/sec)
     &  , fkpm   = 10.*fkph ! background vertical viscosity coefficient (cm**2/sec)

      END MODULE KPP_COM
C****
      MODULE KPPE
!@sum  KPPE contains variables and routines for KPP mixing scheme
!@auth NCAR (modifications by Gavin Schmidt)
c====================== include file "KPP_1D.COM" =====================
c
      USE OCEAN, only : lmo
      USE KPP_COM, only : fkph, fkpm
      USE SW2OCEAN, only : lsrpd,fsr,fsrz,dfsrdz,dfsrdzb
      IMPLICIT NONE
      SAVE
c set max. number of ocean levels from ocean code
!@var km is max number of KPP layers
      INTEGER, PARAMETER :: km=LMO

c     rules for parameter constants
c
c use prefix of "c" for whole real numbers (ie: c57 for 57.0)
c use "m" after prefix to designate negative values (minus sign)
c      (ie: cm7 for -7.0)
c use prefix of "p" for non repeating fractions (ie: p5 for 0.5)
c use prefix of "r" for reciprocals (ie: r3 for 1/3.0)
c combine use of prefix above and "e" for scientific notation, with
c       (ie: c5e4 for 5.0e4, c1em10 for 1.0e-10)
c
      real*8,parameter :: c0=0.0, c1=1.0, c2=2.0, c4=4.0, c5=5.0, c8=8.0
      real*8,parameter :: c16=16.0, c360=360.0
      real*8,parameter :: p25=0.25, p5=0.5, p75=0.75
      real*8,parameter :: epsln=1.0e-20
c
      real*8,parameter :: c24=24.0, c60=60.0, c1440=1440.0
      real*8,parameter :: r24=c1/c24, r60=c1/c60, r1440=c1/c1440
      real*8,parameter :: secday=c1/(c60*c1440)
c
      real*8,parameter :: c1p5 = 1.5 , c3   = 3.     , c35   = 35.
      real*8,parameter :: c10 = 10. , c100 = 100.   , c1000 = 1000.
      real*8,parameter :: c10000 = 10000.
      real*8,parameter :: r3 = c1/c3 , r10  = c1/c10 , r100 = r10/c10
      real*8,parameter :: r1000  = c1/c1000, r10000 = r10*r1000

c====================== include file CVMIX_CLEAN.COM ==================

c     variables used for vertical diffusion
c
c inputs: (set through namelist)
c
!@var fkph   = vertical diffusion coefficient (cm**2/sec) (taken from KPP_COM)
!@var fkpm   = vertical viscosity coefficient (cm**2/sec) (taken from KPP_COM)
!@var bvdc   = background vertical diffusion constant
!@var bvvc   = background vertical viscosity constant
!@var vvclim = vertical viscosity coefficient limit
!@var vdclim = vertical diffusion coefficient limit
c
!@var vvcric = maximum viscosity   due to shear instability    (cm**2/s)
c              (based on local richardson number)
!@var vdcric = maximum diffusivity due to shear instability    (cm**2/s)
c
c arrays used for vertical diffusion in "k-profile" mixing scheme,
c computed in "kmix.F" and "kmixs.F".
c note, that in this scheme temperature diffusvity might be
c different from the diffusivity of all other tracers due to double
c diffusion.
c
      real*8, parameter :: vvcric=50d0, vdcric=50d0
     *     ,vvclim = 1000d0, vdclim = 1000d0

c====================== include file "KMIX_CLEAN.COM" =================
c     Define various parameters and common blocks for kmix vertical-
c     mixing scheme; used in "kmixs.F" subroutines
c
c parameters for several subroutines
c
!@var  epsl    = epsln in "pconst.h" (set in "kmixinit")    = 1.0e-20
!@var  epsilon = nondimensional extent of the surface layer = 0.1
!@var  vonk    = von Karman's constant                      = 0.4
!@var  conc1,conam,concm,conc2,zetam,conas,concs,conc3,zetas
c          = scalar coefficients

      real*8, parameter :: epsl=epsln, epsilon=0.1d0, vonk=0.4d0, conc1
     *     =5d0,conam=1.257d0, concm=8.380d0, conc2=16d0,zetam= -0.2d0
     *     ,conas=-28.86d0,concs=98.96d0,conc3=16d0,zetas= -1d0

c     parameters for subroutine "bldepth"

c to compute depth of boundary layer:
c
!@var  Ricr    = critical bulk Richardson Number            = 0.3 or 1.
!@var  cekman  = coefficient for ekman depth                = 0.7
!@var  cmonob  = coefficient for Monin-Obukhov depth        = 1.0
!@var  concv   = ratio of interior buoyancy frequency to
c               buoyancy frequency at entrainment depth    = 1.8
!@var  hbf     = fraction of bounadry layer depth to
c               which absorbed solar radiation
c               contributes to surface buoyancy forcing    = 1.0
!@var  Vtc     = non-dimensional coefficient for velocity
c               scale of turbulant velocity shear
c               (=function of concv,concs,epsilon,vonk,Ricr)
c
#ifdef OCN_GISS_TURB
      real*8, parameter :: Ricr=1.0d0
#else
      real*8, parameter :: Ricr=0.3d0
#endif
      real*8, parameter :: cekman= 0.7d0, cmonob=1d0, concv
     *     =1.8d0,hbf=1d0
      real*8 Vtc

c     parameters and common arrays for subroutines "kmixinit" and
c     "wscale" to compute turbulent velocity scales:
c
!@var  nni     = number of values for zehat in the look up table
!@var  nnj     = number of values for ustar in the look up table
c
!@var  wmt     = lookup table for wm, the turbulent velocity scale
c               for momentum
!@var  wst     = lookup table for ws, the turbulent velocity scale
c               for scalars
!@var  deltaz  = delta zehat in table
!@var  deltau  = delta ustar in table
!@var  rdeltaz  = recipricol of delta zehat in table
!@var  rdeltau  = recipricol of delta ustar in table
!@var  zmin    = minimum limit for zehat in table (m3/s3)
!@var  zmax    = maximum limit for zehat in table
!@var  umin    = minimum limit for ustar in table (m/s)
!@var  umax    = maximum limit for ustar in table

      integer, parameter :: nni = 890, nnj = 480
c
      real*8 wmt(0:nni+1,0:nnj+1),wst(0:nni+1,0:nnj+1)
      real*8, parameter :: zmin=-4d-7,zmax=0.,umin=0.,umax=4d-2,
     *     deltaz = (zmax-zmin)/(nni+1), deltau = (umax-umin)/(nnj+1),
     *     rdeltaz = 1./deltaz, rdeltau = 1./deltau
c
c     parameters for subroutine "ri_iwmix"
c     to compute vertical mixing coefficients below boundary layer:
c
!@var  Riinfty = local Richardson Number limit
c              for shear instability                      = 0.7
!@var  rRiinfty = 1/Riinfty
c    (note: the vertical mixing coefficients are defined
c    in (m2/s) units in "kmixinit". they are initialized
c    in "kmixbdta" and read via namelist "eddy" in (cm2/s)
c    units using their c-g-s names)
c   (m2/s)                                                 (cm2/s)
!@var difm0   = viscosity max due to shear instability     = vvcric
!@var difs0   = diffusivity ..                             = vdcric
!@var difmiw  = viscosity background due to internal waves = fkpm
!@var difsiw  = diffusivity ..                             = fkph
!@var difmcon = viscosity due to convective instability    = vvclim
!@var difscon = diffusivity ..                             = vdclim
!@var BVSQcon = value of N^2 where the convection coefs first become max
!@var rBVSQcon = 1/BVSQcon
!@var num_v_smooth_Ri = number of vertical smoothings of Ri

      real*8, parameter :: Riinfty= 0.7d0, rRiinfty=1./Riinfty, BVSQcon=
     *     -1d-7,rBVSQcon=1./BVSQcon
      integer, parameter :: num_v_smooth_Ri=1
      real*8, parameter :: difm0   = vvcric * r10000, difs0   = vdcric
     *     * r10000,difmiw  = fkpm   * r10000, difsiw  = fkph   * r10000
     *     ,difmcon = vvclim * r10000, difscon = vdclim * r10000

c     parameters for subroutine "ddmix"
c     to compute additional diffusivity due to double diffusion:
!@var  Rrho0   = limit for double diffusive density ratio
!@var  dsfmax  = maximum diffusivity in case of salt fingering (m2/s)

      real*8, parameter :: Rrho0 = 1.9d0, dsfmax = 0.001d0

c     parameters for subroutine "blmix"
c     to compute mixing within boundary layer:
!@var  cstar   = proportionality coefficient for nonlocal transport
!@var  cg      = non-dimensional coefficient for counter-gradient term

      real*8, parameter :: cstar = 10d0
      real*8 cg

c add variables for depth dependent mixing due to rough topography
!@var diftop diffusion scale for topographic mixing
!@var fz500 diffusion scale as function of height above topography
      real*8, parameter :: diftop=0d0    ! 10d0 * r10000)
      real*8 fz500(km,km)

      END MODULE KPPE


      SUBROUTINE KPPMIX(LDD   , ZE    ,
     $                  zgrid , hwide , kmtj  , Shsq   , dVsq  ,
     $                  ustar , Bo    , Bosol , alphaDT, betaDS,
     $                  dbloc , Ritop , Coriol, byhwide,
     $                  visc  , difs  , dift  , ghats  , hbl
     $                                                 , kbl )

!@sum KPPMIX Main driver subroutine for kpp vertical mixing scheme and
!@+   interface to greater ocean model
c
c     written  by: bill large,    june  6, 1994
c     modified by: jan morzel,    june 30, 1994
c                  bill large,  august 11, 1994
c                  bill large, january 25, 1995 : "dVsq" and 1d code
c     modified for GISS by Gavin Schmidt, march 1998
c              for ModelE                 march 2001
c
      USE KPPE
!@var mdiff number of diffusivities for local arrays
      integer, parameter :: mdiff = 3
c input
      real*8 ZE(0:LMO)      !@var ZE GISS vertical layering (m)
      real*8 zgrid(0:km+1)  !@var zgrid vertical grid (<= 0) (m)
      real*8 hwide(0:km+1)  !@var hwide layer thicknesses    (m)
      real*8 byhwide(0:km+1)!@var byhwide 1/layer thicknesses (1/m)
      integer kmtj    !@var kmtj number of vertical layers on this row
      real*8 Shsq(km) !@var Shsq (local velocity shear)^2  (m/s)^2
      real*8 dVsq(km) !@var dVsq (velocity shear re sfc)^2 (m/s)^2
      real*8 ustar    !@var ustar surface friction velocity (m/s)
      real*8 Bo       !@var Bo surface turbulent buoy. forcing (m^2/s^3)
      real*8 Bosol    !@var Bosol radiative buoyancy forcing (m^2/s^3)
!@var alphaDT alpha * DT  across interfaces (kg/m^3)
      real*8 alphaDT(km)
!@var betaDS beta  * DS  across interfaces (kg/m^3)
      real*8 betaDS(km)
!@var dbloc local delta buoyancy across interfaces (m/s^2)
      real*8 dbloc(km)
!@var Ritop numerator of bulk Richardson Number (m/s)^2
c          Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc
      real*8 Ritop(km)
      real*8 Coriol   !@var Coriol Coriolis parameter            (1/s)
      logical LDD     !@var LDD = TRUE for double diffusion
c output
!@var visc vertical viscosity coefficient (m^2/s)
      real*8 visc(0:km+1)
!@var difs vertical scalar diffusivity (m^2/s)
      real*8 difs(0:km+1)
!@var dift vertical temperature diffusivity (m^2/s)
      real*8 dift(0:km+1)
      real*8 ghats(km)  !@var ghats nonlocal transport (s/m^2)
      real*8 hbl        !@var hbl boundary layer depth (m)
      real*8 byhbl      !@var byhbl 1/boundary layer depth (1/m)
c local
      real*8 bfsfc      !@var bfsfc surface buoyancy forcing (m^2/s^3)
      real*8 ws         !@var ws momentum velocity scale
      real*8 wm         !@var wm scalar   velocity scale
      real*8 caseA      !@var caseA = 1 in case A; =0 in case B
      real*8 stable     !@var stable 1 in stable forcing; 0 in unstable
      real*8 dkm1(mdiff)!@var dkm1 boundary layer difs at kbl-1 level
      real*8 gat1(mdiff)!@var gat1 shape function at sigma=1
!@var dat1 derivative of shape function at sigma=1
      real*8 dat1(mdiff)
      real*8 blmc(km,mdiff)!@var blmc boundary layer mixing coefficients
      real*8 sigma      !@var sigma normalized depth (d / hbl)
      real*8 Rib(2)     !@var Rib bulk Richardson number
      integer kbl       !@var kbl index of first grid level below hbl
      integer kmax  !@var kmax minimum of LSRPD and kmtj, used in swfrac
c local ri_iwmix
      real*8 Rigg       !@var Rigg local richardson number
      real*8 fri,fcon   !@var fri,fcon function of Rig
      real*8 ftop       !@var ftop function of topography
c local enhance
      real*8 delta  !@var delta fraction hbl lies beteen zgrid neighbors
c local wscale
      real*8 zehat      !@var zehat = zeta *  ustar**3
      INTENT (IN) LDD,ZE,zgrid,hwide,kmtj,Shsq,dVsq,ustar,Bo,Bosol
     *            ,alphaDT,betaDS,dbloc,Ritop,Coriol,byhwide
      INTENT (OUT) visc, difs, dift, ghats, hbl, kbl
      INTEGER ki,mr,ka,ku,kl,iz,izp1,ju,jup1,ksave,kn,kt
      REAL*8 ratio,zdiff,zfrac,fzfrac,wam,wbm,was,wbs,u3,bvsq
     *     ,delhat,dvdzup,dvdzdn,viscp,diftp,visch,difsh,f1,bywm,byws
     *     ,sig,a1,a2,a3,gm,gs,gt,dstar,udiff,ufrac,vtsq,r,difsp,difth
     *     ,dkmp5

      kmax = MIN(LSRPD,kmtj)
c compute interior mixing coefficients everywhere, due to constant
c internal wave activity, static instability, and local shear
c instability.
c Zero surface values

c      call ri_iwmix ( km, km+1, kmtj, Shsq, dbloc, zgrid,
c    *                visc, difs , dift  )
c     compute interior viscosity diffusivity coefficients due
c     to shear instability (dependent on a local richardson number),
c     to background internal wave activity, and
c     to static instability (local richardson number < 0).

c     compute interior gradient Ri at all interfaces ki=1,km, (not surfa
ce)
c       use visc(ki=1,km) as temporary storage to be smoothed
c       use dift(ki=1,km) as temporary storage of N^2 = BVSQ
c       set values at bottom and below to nearest value above bottom

      do ki = 1, kmtj
         visc(ki)  = dbloc(ki) * (zgrid(ki)-zgrid(ki+1)) /
     $        ( Shsq(ki) + epsl)
         dift(ki) = dbloc(ki) * byhwide(ki+1)
      end do

c vertically smooth Ri num_v_smooth_Ri times

      do mr = 1,num_v_smooth_Ri
        call z121(visc,kmtj,km)
      end do

      do ki = 1, kmtj

c evaluate f of Vaisala squared    for convection        store in fcon
c evaluate f of   smooth Ri (fri) for shear instability store in fri

         Rigg  = DMAX1( dift(ki) , BVSQcon )
         ratio = DMIN1( (BVSQcon-Rigg) * rBVSQcon , c1 )
         fcon  = (c1 - ratio*ratio)
         fcon  = fcon * fcon * fcon

         Rigg  = DMAX1( visc(ki) , c0 )
         ratio = DMIN1( Rigg * rRiinfty , c1 )
         fri   = (c1 - ratio*ratio)
         fri   = fri  * fri   * fri
c add increased tracer mixing over rough topography. As a first cut,
c assume that topography is 'rough' if kmtj != LMO.
c functional form kv_top = diftop * exp (-z/500) where z is
c distance from the bottom.
c The array FZ500 is zero if kmtj=LMO

          ftop =  FZ500(ki,kmtj)

c evaluate diffusivities and viscosity
c mixing due to internal waves, and shear and static instability

         visc(ki) = (difmiw + fcon * difmcon + fri * difm0)
         difs(ki) = (difsiw + fcon * difscon + fri * difs0
     *                      + ftop * diftop)
         dift(ki) = difs(ki)
      end do

c set surface values to 0.0

      visc(0)    = c0
      dift(0)    = c0
      difs(0)    = c0

c add double diffusion if desired

      if (LDD) call ddmix (alphaDT,betaDS,visc,difs,dift,kmtj)

c Zero values at seafloor and below for blmix

      visc(kmtj:km+1) = c0
      difs(kmtj:km+1) = c0
      dift(kmtj:km+1) = c0

c compute boundary layer mixing coefficients:
c diagnose the new boundary layer depth

c       call bldepth (km, km+1, zgrid, hwide, kmtj, dVsq,
c      $             dbloc, Ritop, ustar, Bo, Bosol, Coriol,
c      $             hbl, bfsfc, stable, caseA, kbl,
c      $             Rib, sigma, wm, ws)

c     the oceanic planetray boundary layer depth, hbl, is determined as
c     the shallowest depth where the bulk richardson number is
c     equal to the critical value, Ricr.
c     Bulk richardson numbers are evaluated by computing velocity and
c     buoyancy differences between values at zgrid(kl) < 0 and surface
c     reference values.
c     in this configuration, the reference values are equal to the
c     values in the surface layer.
c     when using a very fine vertical grid, these values should be
c     computed as the vertical average of velocity and buoyancy from
c     the surface down to epsilon*zgrid(kl).
c     When the bulk richardson number at k exceeds Ricr, hbl is
c     linearly interpolated between grid levels zgrid(k) and zgrid(k-1).
c     The water column and the surface forcing are diagnosed for
c     stable/ustable forcing conditions, and where hbl is relative
c     to grid points (caseA), so that conditional branches can be
c     avoided in later subroutines.

c find bulk Richardson number at every grid level until > Ricr
c
c note: the reference depth is -epsilon/2.*zgrid(k), but the reference
c       u,v,t,s values are simply the surface layer values,
c       and not the averaged values from 0 to 2*ref.depth,
c       which is necessary for very fine grids(top layer < 2m thickness)
c note: max values when Ricr never satisfied are
c       kbl=kmtj and hbl=-zgrid(kmtj)

c indices for array Rib(k), the bulk Richardson number.
      ka = 1
      ku = 2

c initialize hbl and kbl to bottomed out values
      Rib(ka) = 0.0
      kbl    = kmtj
      hbl    = -zgrid(kbl)

      kl = 1
 10   kl = kl + 1

c compute bfsfc = sw fraction at hbf * zgrid
c      call swfrac(zgrid(kl),ZE,kl,kmtj,kmax,bfsfc)  ! inlined

C**** Get actual k (since k is on staggered grid)
      kt = kl + NINT(0.5d0 + SIGN(0.5d0,-(ZE(kl)+zgrid(kl))))

C**** calculate fraction
      if (kt.gt.kmax) then
         bfsfc = 0.
      else
         if (kt.eq.kmax) then
            bfsfc = FSR(kt) + (zgrid(kl)+ZE(kt-1))*dFSRdZB(kt)
         else
            bfsfc = FSR(kt) + (zgrid(kl)+ZE(kt-1))*dFSRdZ(kt)
         end if
      end if

c use caseA as temporary array
         caseA  = -zgrid(kl)

c compute bfsfc= Bo + radiative contribution down to hbf * hbl
         bfsfc  = Bo + Bosol * (1. - bfsfc)

         stable = 0.5 + SIGN( 5d-1, bfsfc )
         sigma  = stable * 1. + (1.-stable) * epsilon

c compute velocity scales at sigma, for hbl= caseA = -zgrid(kl)
c         call wscale(sigma, caseA, ustar, bfsfc,   wm, ws)
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * caseA * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam
         if (ju == nnj .and. wm < wam) wm = wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
         if (ju == nnj .and. ws < was) ws = was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

c compute the turbulent shear contribution to Rib
         bvsq =0.5*
     $        ( dbloc(kl-1) * byhwide(kl) +
     $          dbloc(kl  ) * byhwide(kl+1) )
         Vtsq = - zgrid(kl) * ws * sqrt(abs(bvsq)) * Vtc

c           compute bulk Richardson number at new level, dunder
c           note: Ritop needs to be zero on land and ocean bottom
c           points so that the following if statement gets triggered
c           correctly. otherwise, hbl might get set to (big) negative
c           values, that might exceed the limit for the "exp" function
c           in "swfrac"

         Rib(ku) = Ritop(kl) / (dVsq(kl)+Vtsq+epsl)

c           linearly interpolate to find hbl where Rib = Ricr
         if((kbl.eq.kmtj).and.(Rib(ku).gt.Ricr)) then
            hbl = -zgrid(kl-1) + (zgrid(kl-1)-zgrid(kl)) *
     $           (Ricr - Rib(ka)) / (Rib(ku)-Rib(ka))
            kbl = kl
         else
            ksave = ka
            ka    = ku
            ku    = ksave
            if (kl.lt.kmtj) go to 10
         end if

c find stability and buoyancy forcing for boundary layer
c      call swfrac(-hbl,ZE,kbl-1,kmtj,bfsfc,kmax)  ! inlined

C**** Get actual k (since k is on staggered grid)
      kt = kbl-1+ NINT(0.5d0 + SIGN(0.5d0,-(ZE(kbl-1)-hbl)))

C**** calculate fraction
      if (kt.gt.kmax) then
         bfsfc = 0.
      else
         if (kt.eq.kmax) then
            bfsfc = FSR(kt) + (ZE(kt-1)-hbl)*dFSRdZB(kt)
         else
            bfsfc = FSR(kt) + (ZE(kt-1)-hbl)*dFSRdZ(kt)
         end if
      end if

      bfsfc  = Bo + Bosol * (1. - bfsfc)
      stable = 0.5 + SIGN( 5d-1, bfsfc )
      bfsfc  = bfsfc + stable * epsl !ensures bfsfc never=0

c check hbl limits for hekman or hmonob (NOT USED)

c      if(bfsfc.gt.0.0) then
c         hekman = cekman * ustar / (abs(Coriol)+epsl)
c         hmonob = cmonob * ustar*ustar*ustar
c     $        /(vonk * (bfsfc+epsl) )
c         hlimit = stable     * DMIN1(hekman,hmonob) +
c     $        (stable-1.) * zgrid(km)
c         hbl = DMIN1(hbl,hlimit)
c         hbl = DMAX1(hbl,-zgrid(1))
c      endif
c      kbl = kmtj

c find new kbl
c      do kl=2,kmtj
c         if((kbl.eq.kmtj).and.(-zgrid(kl).gt.hbl)) kbl = kl
c      end do

c find stability and buoyancy forcing for final hbl values
c      call swfrac(-hbl,ZE,kbl-1,kmtj,kmax,bfsfc)

c      bfsfc  = Bo + Bosol * (1. - bfsfc)
c      stable = 0.5 + SIGN( 5d-1, bfsfc )
c      bfsfc  = bfsfc + stable * epsl

c determine caseA and caseB
      caseA  = 0.5 + SIGN(5d-1,-zgrid(kbl) - 0.5*hwide(kbl) - hbl)

      byhbl = 1d0/hbl
c compute boundary layer diffusivities

c       call blmix    (km   , km+1 , mdiff , zgrid, hwide ,
c      $               ustar, bfsfc, hbl  , stable, caseA,
c      $               visc , difs , dift , kbl   ,
c      $               gat1 , dat1 , dkm1 , blmc  , ghats,
c      $               sigma, wm   , ws   )
c     mixing coefficients within boundary layer depend on surface
c     forcing and the magnitude and gradient of interior mixing below
c     the boundary layer ("matching").
c     Caution: if mixing bottoms out at hbl = -zgrid(km) then
c     fictitious layer at km+1 is needed with small but finite width
c     hwide(km+1) (eg. epsl = 1.e-20).

c compute velocity scales at hbl

      sigma = stable * 1.0 + (1.-stable) * epsilon

c      call wscale(sigma, hbl, ustar, bfsfc,   wm, ws)
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * hbl * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam
         if (ju == nnj .and. wm < wam) wm = wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
         if (ju == nnj .and. ws < was) ws = was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

      kn    = int(caseA+epsl) *(kbl -1) +
     $     (1-int(caseA+epsl)) * kbl

c find the interior viscosities and derivatives at hbl
      delhat = 0.5*hwide(kn) - zgrid(kn) - hbl
      R      = 1.0 - delhat * byhwide(kn)
      dvdzup = (visc(kn-1) - visc(kn)) * byhwide(kn)
      dvdzdn = (visc(kn)   - visc(kn+1)) * byhwide(kn+1)
      viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $     R  * (dvdzdn + abs(dvdzdn)) )

      dvdzup = (difs(kn-1) - difs(kn)) * byhwide(kn)
      dvdzdn = (difs(kn)   - difs(kn+1)) * byhwide(kn+1)
      difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $     R  * (dvdzdn + abs(dvdzdn)) )

      dvdzup = (dift(kn-1) - dift(kn)) * byhwide(kn)
      dvdzdn = (dift(kn)   - dift(kn+1)) * byhwide(kn+1)
      diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $     R  * (dvdzdn + abs(dvdzdn)) )

      visch  = visc(kn) + viscp * delhat
      difsh  = difs(kn) + difsp * delhat
      difth  = dift(kn) + diftp * delhat

      f1 = stable * conc1 * bfsfc / (ustar**4+epsl)

      bywm = 1./(wm+epsl)
      byws = 1./(ws+epsl)

      gat1(1) = visch * byhbl * bywm
      dat1(1) = -viscp * bywm + f1 * visch
      dat1(1) = min(dat1(1),0d0)

      gat1(2) = difsh  * byhbl * byws
      dat1(2) = -difsp * byws  + f1 * difsh
      dat1(2) = min(dat1(2),0d0)

      gat1(3) = difth * byhbl * byws
      dat1(3) = -diftp * byws + f1 * difth
      dat1(3) = min(dat1(3),0d0)

      do ki = 1,kbl-1

c compute turbulent velocity scales on the interfaces
         sig     = (-zgrid(ki) + 0.5 * hwide(ki)) * byhbl
         sigma= stable*sig + (1.-stable)*DMIN1(sig,epsilon)

c         call wscale(sigma, hbl, ustar, bfsfc,   wm,  ws)
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * hbl * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam
         if (ju == nnj .and. wm < wam) wm = wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
         if (ju == nnj .and. ws < was) ws = was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

c compute the dimensionless shape functions at the interfaces
         sig = (-zgrid(ki) + 0.5 * hwide(ki)) * byhbl
         a1 = sig - 2.
         a2 = 3.-2.*sig
         a3 = sig - 1.

         Gm = a1 + a2 * gat1(1) + a3 * dat1(1)
         Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
         Gt = a1 + a2 * gat1(3) + a3 * dat1(3)

c compute boundary layer diffusivities at the interfaces
         blmc(ki,1) = hbl * wm * sig * (1. + sig * Gm)
         blmc(ki,2) = hbl * ws * sig * (1. + sig * Gs)
         blmc(ki,3) = hbl * ws * sig * (1. + sig * Gt)

c nonlocal transport term = ghats * <ws>o
         ghats(ki) = (1.-stable) * cg * byws * byhbl
      end do

c find diffusivities at kbl-1 grid level
      sig      =  -zgrid(kbl-1)  * byhbl
      sigma =stable * sig + (1.-stable) * MIN(sig,epsilon)

c      call wscale(sigma, hbl, ustar, bfsfc,   wm, ws)
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * hbl * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam
         if (ju == nnj .and. wm < wam) wm = wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
         if (ju == nnj .and. ws < was) ws = was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

      sig = -zgrid(kbl-1) * byhbl
      a1= sig - 2.
      a2 = 3.-2.*sig
      a3 = sig - 1.
      Gm = a1 + a2 * gat1(1) + a3 * dat1(1)
      Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
      Gt = a1 + a2 * gat1(3) + a3 * dat1(3)
      dkm1(1) = hbl * wm * sig * (1. + sig * Gm)
      dkm1(2) = hbl * ws * sig * (1. + sig * Gs)
      dkm1(3) = hbl * ws * sig * (1. + sig * Gt)

c       call enhance  (km   , km+1 , mdiff , dkm1 , visc  ,
c      $               difs , dift , hbl  , kbl   , zgrid, caseA ,
c      $               blmc ,ghats )
c     enhance the diffusivity at the kbl-.5 interface

      if (kbl.le.kmtj) then
         ki = kbl-1
         delta = (hbl+zgrid(ki)) * byhwide(ki+1)

         dkmp5 = caseA*visc(ki) + (1.-caseA)*blmc(ki,1)
         dstar = (1.-delta)**2 * dkm1(1) + delta**2 * dkmp5
         blmc(ki,1) = (1.-delta) * visc(ki) + delta * dstar

         dkmp5 = caseA*difs(ki) + (1.-caseA)*blmc(ki,2)
         dstar = (1.-delta)**2 * dkm1(2) + delta**2 * dkmp5
         blmc(ki,2) = (1.-delta) * difs(ki) + delta * dstar

         dkmp5 = caseA*dift(ki) + (1.-caseA)*blmc(ki,3)
         dstar = (1.-delta)**2 * dkm1(3) + delta**2 * dkmp5
         blmc(ki,3) = (1.-delta) * dift(ki) + delta * dstar

         ghats(ki) = (1.-caseA) * ghats(ki)
      end if

c combine interior and boundary layer coefficients and nonlocal term

      visc(1:kbl-1)=blmc(1:kbl-1,1)
      difs(1:kbl-1)=blmc(1:kbl-1,2)
      dift(1:kbl-1)=blmc(1:kbl-1,3)
      ghats(kbl:km)=0.

      return
      end subroutine kppmix


      subroutine bldepth(
      ! in:
     &     ze,zgrid,byhwide,lmij,dVsq
     &    ,ustar,Bo,Bosol,dbloc,Ritop
      ! out:
     &    ,rib,hbl,kbl,bf)
!@sum bldepth calculates oceanic planetray boundary layer depth, hbl
!@+     written  by: bill large,    june  6, 1994
!@+     modified by: jan morzel,    june 30, 1994
!@+                  bill large,  august 11, 1994
!@+                  bill large, january 25, 1995 : "dVsq" and 1d code
!@+     modified for GISS by Gavin Schmidt, march 1998
!@+              for ModelE                 march 2001
!@+     modified for GISS by Ye Cheng and Armando Howard, May 2011
 
      USE KPPE, only : lmo,LSRPD,rdeltaz,nni,rdeltau,nnj,wmt,wst,vonk
     &   ,Vtc,Ricr,fsr,dfsrdzb,dfsrdz,zmax,zmin,umin,conc1,epsl,epsilon
      USE DOMAIN_DECOMP_1d, Only: AM_I_ROOT

c input
      real*8 ze(0:lmo)       !@var ze giss vertical layering (m)
      real*8 zgrid(0:lmo+1)  !@var zgrid vertical grid (<= 0) (m)
      real*8 byhwide(0:lmo+1)!@var byhwide 1/layer thicknesses (1/m)
      integer lmij     !@var lmij number of vertical layers on this row
      real*8 dVsq(lmo) !@var dVsq (velocity shear re sfc)^2 (m/s)^2
      real*8 ustar     !@var ustar surface friction velocity (m/s)
      real*8 Bo        !@var Bo surface turbulent buoy. forcing (m^2/s^3)
      real*8 Bosol     !@var Bosol radiative buoyancy forcing (m^2/s^3)
!@var dbloc local delta buoyancy across interfaces (m/s^2)
      real*8 dbloc(lmo)
!@var Ritop numerator of bulk Richardson Number (m/s)^2
      real*8 Ritop(lmo)
      INTENT (IN) ze,zgrid,byhwide,lmij,dVsq,ustar,Bo,Bosol
     *           ,dbloc,Ritop
c output
      real*8 rib(lmo)   !@var rib bulk Richardson number
      real*8 hbl        !@var hbl boundary layer depth (m)
      integer kbl       !@var kbl index of first grid level below hbl
      real*8 bf         !@var bfsfc surface buoyancy forcing (m^2/s^3)
      INTENT (OUT) rib,hbl,kbl,bf
c local
      real*8 bfsfc      !@var bfsfc surface buoyancy forcing (m^2/s^3)
      real*8 rib2(2)    !@var rib2 temperary bulk Richardson number
      real*8 byhbl      !@var byhbl 1/boundary layer depth (1/m)
      real*8 ws         !@var ws momentum velocity scale
      real*8 wm         !@var wm scalar   velocity scale
      real*8 caseA      !@var caseA = 1 in case A; =0 in case B
      real*8 stable     !@var stable 1 in stable forcing; 0 in unstable
      real*8 sigma      !@var sigma normalized depth (d / hbl)
      integer lmax  !@var lmax minimum of LSRPD and lmij, used in swfrac
      real*8 zehat      !@var zehat = zeta *  ustar**3
      INTEGER ki,mr,ka,ku,kl,iz,izp1,ju,jup1,ksave,kt
      REAL*8 zdiff,zfrac,fzfrac,wam,wbm,was,wbs,u3,bvsq
     *     ,delhat,dvdzup,dvdzdn,viscp,diftp,visch,difsh,f1,bywm,byws
     *     ,udiff,ufrac,vtsq

      lmax = MIN(LSRPD,lmij)

c     the oceanic planetray boundary layer depth, hbl, is determined as
c     the shallowest depth where the bulk richardson number is
c     equal to the critical value, Ricr.
c     Bulk richardson numbers are evaluated by computing velocity and
c     buoyancy differences between values at zgrid(kl) < 0 and surface
c     reference values.
c     in this configuration, the reference values are equal to the
c     values in the surface layer.
c     when using a very fine vertical grid, these values should be
c     computed as the vertical average of velocity and buoyancy from
c     the surface down to epsilon*zgrid(kl).
c     When the bulk richardson number at k exceeds Ricr, hbl is
c     linearly interpolated between grid levels zgrid(k) and zgrid(k-1).
c     The water column and the surface forcing are diagnosed for
c     stable/ustable forcing conditions, and where hbl is relative
c     to grid points (caseA), so that conditional branches can be
c     avoided in later subroutines.

c find bulk Richardson number at every grid level until > Ricr
c
c note: the reference depth is -epsilon/2.*zgrid(k), but the reference
c       u,v,t,s values are simply the surface layer values,
c       and not the averaged values from 0 to 2*ref.depth,
c       which is necessary for very fine grids(top layer < 2m thickness)
c note: max values when Ricr never satisfied are
c       kbl=lmij and hbl=-zgrid(lmij)

c indices for array Rib2(k), the temperary bulk Richardson number.
      ka = 1
      ku = 2

c initialize hbl and kbl to bottomed out values
      Rib2(ka) = 0.0
      kbl    = lmij
      hbl    = -zgrid(kbl)

      kl = 1
      rib(1)=0.
 10   kl = kl + 1

c compute bfsfc = sw fraction at hbf * zgrid
c      call swfrac(zgrid(kl),ZE,kl,lmij,lmax,bfsfc)  ! inlined

C**** Get actual k (since k is on staggered grid)
      kt = kl + NINT(0.5d0 + SIGN(0.5d0,-(ZE(kl)+zgrid(kl))))

C**** calculate fraction
      if (kt.gt.lmax) then
         bfsfc = 0.
      else
         if (kt.eq.lmax) then
            bfsfc = FSR(kt) + (zgrid(kl)+ZE(kt-1))*dFSRdZB(kt)
         else
            bfsfc = FSR(kt) + (zgrid(kl)+ZE(kt-1))*dFSRdZ(kt)
         end if
      end if
c     if( AM_I_ROOT() ) then
c        write(301,'(2i4,9e14.4)') kl,kt,1. - bfsfc,Bo,Bosol
c    &   ,Bo + Bosol * (1. - bfsfc)
c     endif 

c use caseA as temporary array
         caseA  = -zgrid(kl)

c compute bfsfc= Bo + radiative contribution down to hbf * hbl
         bfsfc  = Bo + Bosol * (1. - bfsfc)
         bf=bfsfc

         stable = 0.5 + SIGN( 5d-1, bfsfc )
         sigma  = stable * 1. + (1.-stable) * epsilon

c compute velocity scales at sigma, for hbl= caseA = -zgrid(kl)
c         call wscale(sigma, caseA, ustar, bfsfc,   wm, ws)
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * caseA * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam
         if (ju == nnj .and. wm < wam) wm = wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
         if (ju == nnj .and. ws < was) ws = was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

c compute the turbulent shear contribution to Rib
         bvsq =0.5*
     $        ( dbloc(kl-1) * byhwide(kl) +
     $          dbloc(kl  ) * byhwide(kl+1) )
         Vtsq = - zgrid(kl) * ws * sqrt(abs(bvsq)) * Vtc

c           compute bulk Richardson number at new level, dunder
c           note: Ritop needs to be zero on land and ocean bottom
c           points so that the following if statement gets triggered
c           correctly. otherwise, hbl might get set to (big) negative
c           values, that might exceed the limit for the "exp" function
c           in "swfrac"

         Rib2(ku) = Ritop(kl) / (dVsq(kl)+Vtsq+epsl)
         rib(kl)=rib2(ku)

c           linearly interpolate to find hbl where Rib = Ricr
         if((kbl.eq.lmij).and.(Rib2(ku).gt.Ricr)) then
            hbl = -zgrid(kl-1) + (zgrid(kl-1)-zgrid(kl)) *
     $           (Ricr - Rib2(ka)) / (Rib2(ku)-Rib2(ka))
            kbl = kl
         else
            ksave = ka
            ka    = ku
            ku    = ksave
            if (kl.lt.lmij) go to 10
         end if

c find stability and buoyancy forcing for boundary layer
c      call swfrac(-hbl,ZE,kbl-1,lmij,bfsfc,lmax)  ! inlined

C**** Get actual k (since k is on staggered grid)
      kt = kbl-1+ NINT(0.5d0 + SIGN(0.5d0,-(ZE(kbl-1)-hbl)))

C**** calculate fraction
      if (kt.gt.lmax) then
         bfsfc = 0.
      else
         if (kt.eq.lmax) then
            bfsfc = FSR(kt) + (ZE(kt-1)-hbl)*dFSRdZB(kt)
         else
            bfsfc = FSR(kt) + (ZE(kt-1)-hbl)*dFSRdZ(kt)
         end if
      end if

      bfsfc  = Bo + Bosol * (1. - bfsfc)
      stable = 0.5 + SIGN( 5d-1, bfsfc )
      bfsfc  = bfsfc + stable * epsl !ensures bfsfc never=0
      bf=bfsfc

      return
      end subroutine bldepth


      subroutine wscale(sigma, hbl, ustar, bfsfc, wm , ws)

c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c
c     note: the lookup table is only used for unstable conditions
c     (zehat.le.0), in the stable domain wm (=ws) gets computed
c     directly.
      USE KPPE
      INTENT (IN) sigma,hbl,ustar,bfsfc
      INTENT (OUT) wm,ws

c  input
      real*8 sigma      ! normalized depth (d/hbl)
      real*8 hbl        ! boundary layer depth (m)
      real*8 ustar      ! surface friction velocity         (m/s)
      real*8 bfsfc    ! total surface buoyancy flux       (m^2/s^3)
c  output
      real*8 wm,ws ! turbulent velocity scales at sigma
c local
      real*8 zehat           ! = zeta *  ustar**3
      real*8 zdiff,udiff,ufrac,fzfrac,wam,was,wbs,wbm,u3,zfrac
      integer iz,izp1,ju,jup1
c use lookup table for zehat < zmax only; otherwise use
c stable formulae

      zehat = vonk * sigma * hbl * bfsfc

      if (zehat.le.zmax) then
         zdiff  = zehat-zmin
         iz = int( zdiff * rdeltaz )
         iz = min( iz , nni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar-umin
         ju = int( udiff * rdeltau)
         ju = min( ju , nnj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff*rdeltaz - float(iz)
         ufrac = udiff*rdeltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm = (1.-ufrac)* wbm          + ufrac*wam
         if (ju == nnj .and. wm < wam) wm = wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws = (1.-ufrac)* wbs          + ufrac*was
         if (ju == nnj .and. ws < was) ws = was
      else
         u3    = ustar*ustar*ustar
         wm = vonk * ustar * u3 / ( u3 + conc1*zehat )
         ws = wm
      endif

      return
      end subroutine wscale

      subroutine ddmix (alphaDT, betaDS, visc, difs, dift, kmtj)

c     Rrho dependent interior flux parameterization.
c     Add double-diffusion diffusivities to Ri-mix values at blending
c     interface and below.

      USE KPPE
      INTENT (IN) alphaDT,betaDS,kmtj
      INTENT (INOUT) visc, difs, dift
c input
      real*8 alphaDT(km)   ! alpha * DT  across interfaces
      real*8 betaDS(km)    ! beta  * DS  across interfaces
c output
      real*8 visc(0:km+1)  ! interior viscosity           (m^2/s)
      real*8 dift(0:km+1)  ! interior thermal diffusivity (m^2/s)
      real*8 difs(0:km+1)  ! interior scalar  diffusivity (m^2/s)
c local
      real*8 diffdd            ! double diffusion diffusivity scale
      real*8 prandtl           ! prandtl number
      integer kmtj,ki
      real*8 Rrho

      do 100 ki= 1, kmtj

c salt fingering case

         if((alphaDT(ki).gt.betaDS(ki)).and.
     $        (betaDS (ki).gt.0.         )) then

            Rrho       = MIN(alphaDT(ki) / betaDS(ki) , Rrho0)
c           diffdd     = dsfmax*(1.0-((Rrho-1)/(Rrho0-1))**2)**pexp2
            diffdd     =         1.0-((Rrho-1)/(Rrho0-1))**2
            diffdd     = dsfmax*diffdd*diffdd*diffdd
            dift(ki) = dift(ki) + 0.7*diffdd
            difs(ki) = difs(ki) + diffdd

c diffusive convection

         else if ((alphaDT(ki).lt.0.0).and.(betaDS(ki).lt.0.0)
     $           .and.(alphaDT(ki).lt.betaDS(ki)) ) then

            Rrho    = alphaDT(ki) / betaDS(ki)
            diffdd=1.5d-6*9.0*0.101d0*exp(4.6d0*exp(-0.54d0*(1/Rrho-1)))
            prandtl = 0.15*Rrho
            if (Rrho.gt.0.5d0) prandtl = (1.85d0-0.85d0/Rrho)*Rrho
            dift(ki) = dift(ki) + diffdd
            difs(ki) = difs(ki) + prandtl*diffdd
         endif
 100  continue
      return
      end subroutine ddmix

      subroutine kmixinit(ZE)
c     initialize some constants for kmix subroutines, and initialize
c     for kmix subroutine "wscale" the 2D-lookup table for wm and ws
c     as functions of ustar and zetahat (=vonk*sigma*hbl*bfsfc).

      USE KPPE
c local
      real*8 zehat                        ! = zeta *  ustar**3
      real*8 zeta                       ! = stability parameter d/L
      real*8 ZE(0:LMO)                  ! GCM vertical grid
      integer lbot,j,i,l
      real*8 usta

c define some non-dimensional constants and
c the vertical mixing coefficients in m-k-s units

c      epsl    = epsln
      Vtc     = concv * sqrt(0.2/concs/epsilon) / vonk**2 / Ricr
      cg      = cstar * vonk * (concs * vonk * epsilon)**(1./3.)

c      difm0   = vvcric * r10000
c      difs0   = vdcric * r10000
c      difmiw  = fkpm   * r10000
c      difsiw  = fkph   * r10000
c      difmcon = vvclim * r10000
c      difscon = vdclim * r10000

c construct the wm and ws lookup tables

c      deltaz = (zmax-zmin)/(nni+1)
c      deltau = (umax-umin)/(nnj+1)
c      rdeltaz = 1./deltaz
c      rdeltau = 1./deltau

      do 100 i=0,nni+1
         zehat = deltaz*(i) + zmin
         do 90 j=0,nnj+1
            usta = deltau*(j) + umin
            zeta = zehat/(usta**3+epsl)

            if(zehat.ge.0.) then
               wmt(i,j) = vonk*usta/(1.+conc1*zeta)
               wst(i,j) = wmt(i,j)
            else
               if(zeta.gt.zetam) then
                  wmt(i,j) = vonk* usta * sqrt(sqrt(1.-conc2*zeta))
               else
                  wmt(i,j) = vonk* (conam*usta**3-concm*zehat)**(1./3.)
               endif
               if(zeta.gt.zetas) then
                  wst(i,j) = vonk* usta * sqrt(1.-conc3*zeta)
               else
                  wst(i,j) = vonk* (conas*usta**3-concs*zehat)**(1./3.)
               endif
            endif
 90      continue
 100  continue

c  set up depth dependent factor for mixing due to rough topography
c  assume FZ500 = exp (-z/500)
c  for two levels above bottom, starting at level km/2

      FZ500=0.
      do lbot=LMO/2,LMO-1
         do l=lbot-2,lbot-1
            FZ500(l,lbot)=exp(-(ZE(lbot)-ZE(l))/500.)
         end do
      end do

      return
      end subroutine kmixinit

      subroutine swfrac(z,ZE,k,kmtj,kmax,bfsfc)
!@sum swfrac Calculate fraction of solar energy penetrating to depth z
!@+   using linear interpolation for depths between grid levels
!@+   There is a slight error since 'z' is scaled by free surface height
!@+   and ZE is not.
C****
!@var  k is grid box containing point
!@var  kmtj is number of grid points in column
!@var  kmax is the number of grid points that recieve solar radiation
C****
      USE KPPE
      INTENT (IN) z,k,kmtj,kmax,ZE
      INTENT (OUT) bfsfc
      real*8 z,bfsfc
      real*8 ZE(0:LMO)
      integer k,l,kmtj,kmax,kt

C**** Get actual k (since k is on staggered grid)
      kt = k + NINT(0.5d0 + SIGN(0.5d0,-(ZE(k)+z)))

C**** calculate fraction
      if (kt.gt.kmax) then
         bfsfc = 0.
      else
         if (kt.eq.kmax) then
            bfsfc = FSR(kt) + (z+ZE(kt-1))*dFSRdZB(kt)
         else
            bfsfc = FSR(kt) + (z+ZE(kt-1))*dFSRdZ(kt)
         end if
      end if

      return
      end subroutine swfrac

      subroutine z121 (v,kmtj,km)
!@sum z121 Apply 121 smoothing in k to 2-d array V(k=1,km)
!@+   top (0) value is used as a dummy
!@+   bottom (km+1) value is set to input value from above.
      IMPLICIT NONE
      REAL*8, PARAMETER :: p5=5d-1, p25=2.5d-1
      INTEGER, INTENT (IN) :: kmtj,km
      REAL*8, INTENT (INOUT) :: V(0:km+1)  ! 2-D array to be smoothed
      INTEGER K
      REAL*8 tmp

      V(0)      =  p25 * V(1)
      V(kmtj+1) =        V(kmtj)

      do k=1,kmtj
         tmp      =  V(k)
         V(k)   =  V(0)  + p5 * V(k) + p25 * V(k+1)
         V(0)   =  p25 * tmp
      end do
      return
      end subroutine z121

C**** Is this still necessary now that fluxes are saved?
      SUBROUTINE KVINIT
!@sum KVINIT Initialise KMIX and save pre-source term surface values of
!@+   enthalpy, mass and horizontal gradients for kppmix calculation
      USE OCEAN, only : im,jm,lmo,gxmo,sxmo,gymo,symo,g0m,s0m,mo,uo,vo
     &     ,uod,vod
#ifdef TRACERS_OCEAN
     *     ,trmo,txmo,tymo
#endif
      USE KPP_COM, only : G0M1,MO1,GXM1,GYM1,S0M1,SXM1,SYM1,UO1,VO1
     *     ,UOD1,VOD1,lsrpd
#ifdef TRACERS_OCEAN
     *     ,trmo1,txmo1,tymo1
#endif
!      use domain_decomp_1d, only : grid, get
      use domain_decomp_1d, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      INTEGER I,J
      integer :: j_0,j_1

      call getDomainBounds(grid, j_strt=j_0, j_stop=j_1)
C**** Save surface values
      DO J=J_0,J_1
         DO I=1,IM
            G0M1(I,J,1:LSRPD) = G0M(I,J,1:LSRPD)
            GXM1(I,J) = GXMO(I,J,1)
            GYM1(I,J) = GYMO(I,J,1)
            S0M1(I,J) = S0M(I,J,1)
            SXM1(I,J) = SXMO(I,J,1)
            SYM1(I,J) = SYMO(I,J,1)
            MO1(I,J)  = MO(I,J,1)
            UO1(I,J)  = UO(I,J,1)
            VO1(I,J)  = VO(I,J,1)
            UOD1(I,J)  = UOD(I,J,1)
            VOD1(I,J)  = VOD(I,J,1)
#ifdef TRACERS_OCEAN
            TRMO1(:,I,J) = TRMO(I,J,1,:)
            TXMO1(:,I,J) = TXMO(I,J,1,:)
            TYMO1(:,I,J) = TYMO(I,J,1,:)
#endif
         END DO
      END DO
      RETURN
      END SUBROUTINE KVINIT

      SUBROUTINE OCONV
C****
!@sum  OCONV does vertical mixing using coefficients from KPP scheme
!@auth Gavin Schmidt/Gary Russell
!@ver  2009/08/25
C****
      USE CONSTANT, only : grav,omega,UNDEF_VAL
      USE OCEAN, only : im,jm,lmo,g0m,s0m,gxmo,sxmo,symo,gymo,szmo,gzmo
     *     ,ogeoz,hocean,ze,bydxypo,mo,sinpo,dts,lmm,lmv,lmu,ramvs
     *     ,dxypo,cosic,sinic,uo,vo,uod,vod,ramvn,bydts, IVNP,kpl
      USE ODIAG, only : oijl=>oijl_loc,oij=>oij_loc,oijmm
     *     ,ij_hbl,ij_hblmax,ij_bo,ij_bosol,ij_ustar,ijl_kvm,ijl_kvg
     *     ,ijl_kvx,ijl_wgfl,ijl_wsfl,ol,l_rho,l_temp,l_salt  !ij_ogeoz
     *     ,ij_mld,ij_mldmax
#ifdef OCN_GISS_TURB
     *     ,ijl_ri,ijl_rrho,ijl_bv2,ijl_otke,ijl_kvs,ijl_kvc,ijl_buoy
#endif
#ifdef OCN_GISS_SM
     *     ,ijl_fvb
#endif
      USE KPP_COM, only : g0m1,s0m1,mo1,gxm1,gym1,sxm1,sym1,uo1,vo1
     &     ,uod1,vod1
      USE KPP_COM, only : use_tdiss,tdiss,tdiss_n
#ifdef OCN_GISS_TURB
      USE GISSMIX_COM, only : otke,rhobot,exya,ut2a,taubx,tauby
#endif
      USE OFLUXES, only : oRSI, oSOLARw,oSOLARi, oDMUA,oDMVA,oDMUI,oDMVI
      USE SW2OCEAN, only : fsr,lsrpd
      Use DOMAIN_DECOMP_1d, Only: GETDomainBounds, HALO_UPDATE, NORTH,
     *    SOUTH, HALO_UPDATE_BLOCK, AM_I_ROOT, GLOBALSUM
      USE OCEANR_DIM, only : grid=>ogrid
      USE OCEAN, ONLY : GXXMO,GYYMO,GZZMO,GXYMO,SXXMO,SYYMO,SZZMO,SXYMO
      USE OCEAN, ONLY : USE_QUS,NBYZM,I1YZM,I2YZM,DZO
#ifdef ENHANCED_DEEP_MIXING
      USE KPP_COM, only : fkph
      USE CONSTANT, only : pi
#endif
#ifdef TRACERS_OCEAN
      use ocean, only : ntrtrans,motr
      USE OCEAN, only : trmo,txmo,tymo,tzmo
      USE OCEAN, only : txxmo,tyymo,tzzmo,txymo
      Use KPP_COM, Only: trmo1,txmo1,tymo1
      Use ODIAG, Only: toijl=>toijl_loc,toijl_wtfl
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
#endif
#ifdef OCN_GISS_SM
      use giss_sm_com, only : au_sm,av_sm,rx_sm,ry_sm
     &                       ,gx_sm,gy_sm,sx_sm,sy_sm,fvb_sm
     &                       ,p3d,rho3d
     &                       ,rx,ry,gx,gy,sx,sy
#endif

      IMPLICIT NONE

      LOGICAL*4 QPOLE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     &     UT,VT, UTD,VTD
      Real*8 MML0(LMO), MML(LMO), UL0(LMO,IM+2),UL(LMO,IM+2),
     &      ULD0(LMO,IM+2),ULD(LMO,IM+2),
     *      G0ML0(LMO),G0ML(LMO), GXML(LMO),GYML(LMO),GZML(LMO),
     *      S0ML0(LMO),S0ML(LMO), SXML(LMO),SYML(LMO),SZML(LMO),
     *      GXXML(LMO),GYYML(LMO),GXYML(LMO),GZZML(LMO),
     *      SXXML(LMO),SYYML(LMO),SXYML(LMO),SZZML(LMO),
     *      BYMML(LMO),DTBYDZ(LMO),BYDZ2(LMO),RAVM(IM+2),RAMV(IM+2),
     *      BYMML0(LMO),MMLT(LMO),BYMMLT(LMO),
     *      AKVM(0:LMO+1),AKVG(0:LMO+1),AKVS(0:LMO+1),GHATM(LMO,IM+2),
     *      GHATG(LMO),GHATS(LMO),FLG(LMO),FLS(LMO),TXY,
     *      FLDUM(LMO),GHATDUM(LMO),AKVC(0:LMO+1),GHATC(LMO),
     *      DTP4UV(LMO,IM+2),DTP4G(LMO),DTP4S(LMO)
#ifdef OCN_GISS_SM
      real*8, dimension(lmo) :: rxl,ryl,gxl,gyl,sxl,syl
#endif
      INTEGER LMUV(IM+2)
C**** CONV parameters: BETA controls degree of convection (default 0.5).
      REAL*8, PARAMETER :: BETA=5d-1, BYBETA=1d0/BETA
C**** KPP variables
      REAL*8, PARAMETER :: epsln=1d-20
      REAL*8 zgrid(0:LMO+1),hwide(0:LMO+1),Shsq(LMO),dVsq(LMO)
     *     ,talpha(LMO),sbeta(LMO),galpha(LMO)
     &     ,dbloc(LMO),dbsfc(LMO),Ritop(LMO)
     *     ,alphaDT(LMO),betaDS(LMO),alphaDG(LMO)
     &     ,ghat(LMO),byhwide(0:LMO+1)
      REAL*8 G(LMO),S(LMO),TO(LMO),BYRHO(LMO),RHO(LMO),PO(LMO)
      REAL*8, DIMENSION(0:LMO) :: POE
      REAL*8 UKJM(LMO,IM+2)  !  ,UKM(LMO,4,IM,2:JM-1),OLJ(3,LMO,JM)
      REAL*8, DIMENSION(LMO,4,IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     UKM,UKMD
      REAL*8
     *     OLJ(3,LMO,grid%J_STRT_HALO:grid%J_STOP_HALO),OLtemp(LMO,3)
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     &     MA,KLEN,GSAVE3D,SSAVE3D
      REAL*8, DIMENSION(0:LMO,IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     AKVG3D,AKVS3D,AKVC3D,FLG3D,FLS3D
#ifdef OCN_GISS_SM
     &    ,DTP4G3D,DTP4S3D
#endif
      REAL*8, DIMENSION(LMO,IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     DZ3D
      REAL*8 bydz,rhoe
#ifdef TRACERS_OCEAN
      REAL*8 FLT3D(0:LMO,tracerlist%getsize(),IM,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8 TRSAVE3D(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO,
     &                  tracerlist%getsize())
#endif

#ifdef OCN_GISS_TURB
      LOGICAL, PARAMETER :: LDD = .true.
#else
      LOGICAL, PARAMETER :: LDD = .false.
#endif
      INTEGER I,J,K,L,LMIJ,KMUV,IM1,ITER,NSIGG,NSIGS,KBL,II,N,IB
      REAL*8 CORIOL,UISTR,VISTR,U2rho,DELTAM,DELTAE,DELTASR,ANSTR
     *     ,ZSCALE,HBL,HBLP,Ustar,BYSHC,B0,Bosol,R,R2,DTBYDZ2,DM
     *     ,RHOM,RHO1,Bo,DELTAS
      REAL*8 VOLGSP,ALPHAGSP,BETAGSP,TEMGSP,SHCGS,TEMGS,VOLGS
#ifdef ENHANCED_DEEP_MIXING
!@var kvextra an array prescribing background diffusivity as a function of depth
      REAL*8, DIMENSION(LMO) :: KVEXTRA
#endif
!@var kvtdiss profile of diffusivity (m2/s) output by get_kvtdiss
!@var m1d layer mass (kg/m2) input to get_kvtdiss
      real*8, dimension(lmo) :: kvtdiss,m1d

      integer :: j_0,j_1,j_0s,j_1s,j_0h
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      logical ::
     &     adjust_zslope_using_flux,
     &     relax_subgrid_zprofile,
     &     extra_slope_limitations,
     &     mix_tripled_resolution
#ifdef OCN_GISS_TURB
      real*8, parameter :: cd=3.d-3 !@var cd dry drag coeff.
      REAL*8 bf          !@var bf surface buoyancy forcing
c     REAL*8 omfrac      !@var omfrac 1 - fraction of Bosol penetrated
      REAL*8 u2b         !@var u2b velocity squared at zgrid(lmij)
      real*8 ut2         !@var ut2 unresolved bottom velocity squared (m/s)^2
      real*8 exy         !@var exy tidal power input used to tidal diffusivities
      real*8 ustarb2     !@var ustarb2 velocity squared at zg(n)

      REAL*8 uob         !@var uob bottom u-component of velocity at zgrid(lmij) (m/s)
      REAL*8 vob         !@var vob bottom u-component of velocity at zgrid(lmij) (m/s)
      REAL*8 uodb        !@var uob at v velocity point
      REAL*8 vodb        !@var vob at u velocity point
      REAL*8 u2bx        !@var u2b at u velocity point
      REAL*8 u2by        !@var u2b at v velocity point
      real*8 ut2x        !@var ut2x ut2 at the u velocity point
      real*8 ut2y        !@var ut2y ut2 at the v velocity point
c     real*8 unp,vnp
c     integer ier

      REAL*8 rib(lmo)    !@var rib bulk Richardson number
c     REAL*8 wtnl(lmo)   !@var wtnl non-local term of wt
c     REAL*8 wsnl(lmo)   !@var wsnl non-local term of ws
c     real*8 ga          !@var ga grav*alpha
c     real*8 gb          !@var gb grav*beta
c     real*8 buoynl(lmo) !@var buoynl non-local part of buoyancy flux (m^2/s^3)
      REAL*8 ri(0:lmo+1)   !@var ri Rchardson number
      REAL*8 rrho(0:lmo+1) !@var rrho salt to head density ratio
      real*8 bv2(0:lmo+1)  !@var bv2 Brunt Vaisala frequency squared (1/s**2)
      real*8 buoy(0:lmo+1) !@var buoy buoyancy flux (m**2/s**3)
      REAL*8 e(lmo)        !@var e ocean turbulent kinetic energy (m/s)**2
      integer strait
#endif
#ifdef OCN_GISS_SM
c     REAL*8 dfvgdz(lmo) !@var d(fvg)/dz, fvg the counter-gradient part of wg
c     REAL*8 dfvsdz(lmo) !@var d(fvs)/dz, fvs the counter-gradient part of ws
      REAL*8 fvb(0:lmo+1)    !@var fvb the counter-gradient part of wb
      REAL*8 fvg(0:lmo+1)    !@var fvg the counter-gradient part of wg
      REAL*8 fvs(0:lmo+1)    !@var fvs the counter-gradient part of ws
      REAL*8 p4uv(lmo,2) !@var p4uv coeff P4 in the U,V eqn
      REAL*8 uc(lmo)     !@var x-component of velocity at cell center
      REAL*8 vc(lmo)     !@var y-component of velocity at cell center
      REAL*8, PARAMETER :: wta=exp(-1.d0/240.d0) ! average over 5 days
c     REAL*8, PARAMETER :: wta1=exp(-1.d0/1440.d0) ! average over 30 days
      REAL*8 g1,s1,p1,bydz
#endif
#ifdef TRACERS_OCEAN
      Real*8 TRML(LMO,tracerlist%getsize()),TRML1(tracerlist%getsize()),
     &   TZML(LMO,tracerlist%getsize()),TZZML(LMO,tracerlist%getsize()),
     &   DELTATR(tracerlist%getsize()),GHATT(LMO,tracerlist%getsize()),
     *   FLT(LMO,tracerlist%getsize()),DTP4TR(LMO,tracerlist%getsize())
      type(ocn_tracer_entry), pointer :: entry
      REAL*8, DIMENSION(LMO) :: TXML,TYML,TXXML,TYYML,TXYML,
     &     DTBYDZ_TR,BYDZ2_TR
      INTEGER NSIGT
      REAL*8 :: DFLUX,MINRAT ! for GHATT limits
#endif
      real*8 ptdd,ptdm,ptd(lmo),mld
      integer kmld,ip1,lmix

      call getDomainBounds(grid, j_strt=j_0, j_stop=j_1,
     &                j_strt_skp=j_0s, j_stop_skp=j_1s,
     * HAVE_SOUTH_POLE=HAVE_SOUTH_POLE, HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      call getDomainBounds(grid,j_strt_halo=j_0h)

      DTP4G=0.d0; DTP4S=0.d0; DTP4UV=0.d0;
#ifdef OCN_GISS_SM
      DTP4G3D=0.d0; DTP4S3D=0.d0
#endif


C**** initialise diagnostics saved over quarter boxes and longitude
      OLJ = 0.
      OLtemp=0.d0
C**** Load UO,VO into UT,VT.  UO,VO will be updated, while UT,VT
C**** will be fixed during convection.
      call halo_update (grid, VO, from=south)
      call halo_update (grid, UOD, from=south)
      ukm = 0
      ukmd = 0
      DO L=1,LMO
        UT(:,:,L) = UO(:,:,L)
        VT(:,:,L) = VO(:,:,L)
        UTD(:,:,L) = UOD(:,:,L)
        VTD(:,:,L) = VOD(:,:,L)
      END DO

      kvtdiss = 0.

      if(use_qus==1) then
        adjust_zslope_using_flux = .false.
        relax_subgrid_zprofile = .true.
        extra_slope_limitations = .false. ! not needed in QUS flux limiter
        ! At some point, the full kpp scheme could operate using
        ! the tripled vertical resolution.  For now, this option
        ! only applies to the adjustment of the vertical moments.
        mix_tripled_resolution = .not. relax_subgrid_zprofile
      else
        adjust_zslope_using_flux = .true.
        relax_subgrid_zprofile = .false.
        extra_slope_limitations = .true.
        mix_tripled_resolution = .false.
      endif

      if(mix_tripled_resolution) then
        ! save the state that is contemporaneous with the subgrid profiles
        gsave3d = g0m
        ssave3d = s0m
#ifdef TRACERS_OCEAN
        trsave3d = trmo
#endif
      endif

      if(extra_slope_limitations) then
        ! The vertical slopes are adjusted after the mixing and
        ! are not used during mixing, so this section exists
        ! merely to reproduce previous results.
        do l=1,lmo
        do j=j_0s,j_1s
        do ib=1,nbyzm(j,l)
        do i=i1yzm(ib,j,l),i2yzm(ib,j,l)
          if ( abs(SZMO(I,J,L)) > S0M(I,J,L) )
     *         SZMO(I,J,L) = sign(S0M(I,J,L),SZMO(I,J,L)+0d0)
#ifdef TRACERS_OCEAN
          DO N = 1,tracerlist%getsize()
            entry=>tracerlist%at(n)
            if (entry%t_qlimit) then
              if ( abs(TZMO(I,J,L,N)) > TRMO(I,J,L,N) )
     *             TZMO(I,J,L,N) = sign(TRMO(I,J,L,N),TZMO(I,J,L,N))
            end if
C****
          END DO
#endif
        enddo
        enddo
        enddo
        enddo
      endif ! extra_slope_limitations

#ifdef ENHANCED_DEEP_MIXING
      do l=1,lmo-1
        ! option 1: B-L profile
        kvextra(l) = .75d-5 +
     &     (1.325d-4-.75d-5)*(atan((ze(l)-2000d0)/200d0)/pi + 0.5d0)
        ! all options: subtract default background value
        kvextra(l) = kvextra(l) - fkph*1d-4
#ifdef OCN_GISS_TURB
        Untested.  Subtract something other than fkph.
#endif
      enddo
#endif

#ifdef OCN_GISS_SM
      call getdomainbounds(grid,
     &     j_strt=j_0, j_stop=j_1,
     &     j_strt_skp=j_0s, j_stop_skp=j_1s,
     &     have_north_pole=have_north_pole)

      p3d=0.d0
      do l=1,lmo
      do j=j_0,j_1
      do n=1,nbyzm(j,l)
      do i=i1yzm(n,j,l),i2yzm(n,j,l)
        if(l.eq.1) then
          p3d(i,j,l)=5d-1*MO(I,J,1)*GRAV
        else
          p3d(i,j,l)=p3d(i,j,l-1)+5d-1*GRAV*(MO(I,J,L-1)+MO(I,J,L))
        endif
      enddo
      enddo
      enddo
      enddo
      if(have_north_pole) then
        do l=1,lmo
          p3d(2:im,jm,l) = p3d(1,jm,l)
        enddo
      endif

      ! find 3-d rho
      rho3d=0.d0
      do l=1,lmo
      do j=j_0,j_1
      do n=1,nbyzm(j,l)
      do i=i1yzm(n,j,l),i2yzm(n,j,l)
        g1=G0M(I,J,l)/(MO(I,J,l)*DXYPO(J))
        s1=S0M(I,J,l)/(MO(I,J,l)*DXYPO(J))
        p1=p3d(i,j,l)
        rho3d(i,j,l)= 1.d0/volgsp(g1,s1,p1)
      enddo
      enddo
      enddo
      enddo
      if(have_north_pole) then
        do l=1,lmo
          rho3d(2:im,jm,l) = rho3d(1,jm,l)
        enddo
      endif
      call get_gradients0(mo,rho3d,0,rx,ry)
      call get_gradients0(mo,g0m,1,gx,gy)
      call get_gradients0(mo,s0m,1,sx,sy)
#endif

C****
C**** Outside loop over J
C**** Processes are checked and applied on every horizontal quarter box.
C****
      call halo_update (grid,   VO1, from=south)
      call halo_update (grid, oDMVI, from=south)
#ifdef OCN_GISS_TURB
      call halo_update (grid,  ut2a, from=north)
      taubx=0.d0; tauby=0.d0
#endif
      DO 790 J=j_0s,j_1
C**** coriolis parameter, defined at tracer point
      Coriol = 2d0*OMEGA*SINPO(J)

      QPOLE = J.EQ.JM
      IF (QPOLE)  THEN
C****
C**** Define parameters and currents at the North Pole
C****
      LMIJ=LMM(1,JM)
      KMUV=IM+2
      DO I=1,IM
        LMUV(I) = LMV(I,JM-1)
        RAVM(I) = 1d0/IM
        RAMV(I) = RAMVS(JM)
        UL0(1,I) = VT(I,JM-1,1)
        UL(1,I)  = VO1(I,JM-1)
        DO L=2,LMIJ
          UL0(L,I) = VT(I,JM-1,L)
          UL(L,I)  = UL0(L,I)
        END DO
      END DO
      LMUV(IM+1) = LMM(1,JM)
      RAVM(IM+1) = 1d0
      RAMV(IM+1) = 1d0
      UL0(1,IM+1) = UT(IM,JM,1)
      UL(1,IM+1)  = UO1(IM,JM)
      LMUV(IM+2)  = LMM(1,JM)
      RAVM(IM+2)  = 1d0
      RAMV(IM+2)  = 1d0
      UL0(1,IM+2) = UT(IVNP,JM,1)
      UL(1,IM+2)  = UO1(IVNP,JM)
      MML0(1)     = MO(1,JM,1)*DXYPO(JM)
      BYMML0(1)   = 1d0/MML0(1)
      DTBYDZ(1)   = DTS/MO(1,JM,1)
      G0ML0(1)    = G0M(1,JM,1)
      S0ML0(1)    = S0M(1,JM,1)
      MMLT(1)     = MO1(1,JM)*DXYPO(JM)
      BYMMLT(1)   = 1d0/MMLT(1)
      S0ML(1)     = S0M1(1,JM)
      DO L=2,LMIJ
        UL0(L,IM+1) = UT(IM,JM,L)
        UL(L,IM+1)  = UL0(L,IM+1)
        UL0(L,IM+2) = UT(IVNP,JM,L)
        UL(L,IM+2)  = UL0(L,IM+2)
        MML0(L)     = MO(1,JM,L)*DXYPO(JM)
        BYMML0(L)   = 1d0/MML0(L)
        DTBYDZ(L)   = DTS/MO(1,JM,L)
        BYDZ2(L-1)  = 2d0/(MO(1,JM,L-1)+MO(1,JM,L))
        G0ML0(L)    = G0M(1,JM,L)
        S0ML0(L)    = S0M(1,JM,L)
        MMLT(L)     = MML0(L)
        BYMMLT(L)   = BYMML0(L)
        S0ML(L)     = S0ML0(L)
      END DO
      DO L=1,MIN(LSRPD,LMIJ)
        G0ML(L) = G0M1(1,JM,L)
      END DO
      DO L=LSRPD+1,LMIJ
        G0ML(L) = G0ML0(L)
      END DO
#ifdef TRACERS_OCEAN
      TRML1(:)=TRMO1(:,1,JM)
      DO L = 1,LMIJ
        TRML(L,:) = TRMO(1,JM,L,:)
      END DO
#endif
C**** Calculate whole box turbulent stresses sqrt(|tau_0|/rho_0)
C**** except for RHO factor and sqrt: U2rho = Ustar**2 * rho
C**** DM[UV]A are defined on t-grid. DM[UV]I defined on u,v grid
C**** Note that rotational ice stresses are not calculated at the pole
C**** DMUA/I now defined over whole box, not just surface type
      UISTR = 0
      VISTR = 0
      DO I=1,IM
        UISTR = UISTR + oDMVI(I,JM-1)*COSIC(I)
        VISTR = VISTR - oDMVI(I,JM-1)*SINIC(I)
      END DO
      UISTR = UISTR/IM
      VISTR = VISTR/IM
      U2rho = SQRT(
     *     (oDMUA(1,JM) + UISTR)**2+
     *     (oDMVA(1,JM) + VISTR)**2)*BYDTS
C**** Calculate surface mass, salt and heat fluxes
      DELTAM = (MO(1,JM,1) -  MO1(1,JM))*BYDTS
      DELTAE = (G0ML0(1) - G0ML(1))*BYDXYPO(JM)*BYDTS
      DELTAS = (S0ML0(1) - S0ML(1))*BYDXYPO(JM)*BYDTS
      DELTASR= (oSOLARw(1,JM)*(1d0-oRSI(1,JM))+oSOLARi(1,JM)
     *     *oRSI(1,JM))*BYDTS               ! W/m^2
#ifdef TRACERS_OCEAN
      DELTATR(:) = (TRML(1,:)-TRML1(:))*BYDXYPO(JM)*BYDTS  !  kg/m2*s
#endif
      KPL(1,JM) = 1  ! Initialize mixed layer depth
      I=1
      GOTO 500
      END IF
C****
C**** Define parameters and currents away from the poles
C****
      IM1=IM
      I=1
  210 IF(LMM(I,J).LE.1)  GO TO 730
      LMIJ=LMM(I,J)
      KMUV=4
      LMUV(1) = LMU(IM1,J)
      LMUV(2) = LMU(I  ,J)
      LMUV(3) = LMV(I,J-1)
      LMUV(4) = LMV(I,J  )
      MML0(1)   = MO(I,J,1)*DXYPO(J)
      BYMML0(1) = 1 / MML0(1)
      DTBYDZ(1) = DTS/MO(I,J,1)
      MMLT(1)   = MO1(I,J)*DXYPO(J)
      BYMMLT(1) = 1 / MMLT(1)
      DO L=2,LMIJ
        MML0(L)   = MO(I,J,L)*DXYPO(J)
        BYMML0(L) = 1d0/MML0(L)
        DTBYDZ(L) = DTS/MO(I,J,L)
        BYDZ2(L-1)= 2d0/(MO(I,J,L)+MO(I,J,L-1))
        MMLT(L)   = MML0(L)
        BYMMLT(L) = BYMML0(L)
      END DO
C**** Calculate whole box turbulent stresses sqrt(|tau_0|/rho_0)
C**** except for RHO factor and sqrt: U2rho = Ustar**2 * rho
C**** DM[UV]A are defined on t-grid. DM[UV]I defined on u,v grid
C**** DMUA/I now defined over whole box, not just surface type
      UISTR = 0
      VISTR = 0
      IF (oRSI(I,J).gt.0) THEN
        IF (LMU(I,J).gt.0)   UISTR = UISTR + oDMUI(I,J)
        IF (LMU(IM1,J).gt.0) UISTR = UISTR + oDMUI(IM1,J)
        ANSTR=1.-SIGN(0.25,LMU(I,J)-0.5)-SIGN(0.25,LMU(IM1,J)-0.5)
        UISTR = UISTR*ANSTR
        IF (LMV(I,J).gt.0)   VISTR = VISTR + oDMVI(I,J)
        IF (LMV(I,J-1).gt.0) VISTR = VISTR + oDMVI(I,J-1)
        ANSTR=1.-SIGN(0.25,LMV(I,J)-0.5)-SIGN(0.25,LMV(I,J-1)-0.5)
        VISTR = VISTR*ANSTR
      END IF
      U2rho = SQRT((oDMUA(I,J) + UISTR)**2 +
     *             (oDMVA(I,J) + VISTR)**2)*BYDTS
C**** Calculate surface mass flux and Solar forcing
      DELTAM = (MO(I,J,1) -  MO1(I,J))*BYDTS ! kg/m^2 s
      DELTASR = (oSOLARw(I,J)*(1d0-oRSI(I,J))+oSOLARi(I,J)*oRSI(I,J))
     *     *BYDTS               ! W/m^2
      KPL(I,J) = 1  ! Initialize mixed layer depth

      RAVM(1:4) = .5
      RAMV(1:2) = .5
      RAMV(3) = RAMVS(J)
      RAMV(4) = RAMVN(J)
      UL0(1,1) = UT(IM1,J,1)
      UL0(1,2) = UT(I  ,J,1)
      UL0(1,3) = VT(I,J-1,1)
      UL0(1,4) = VT(I,J  ,1)

      ULD0(1,1) = VTD(IM1,J,1)
      ULD0(1,2) = VTD(I  ,J,1)
      ULD0(1,3) = UTD(I,J-1,1)
      ULD0(1,4) = UTD(I,J  ,1)

      G0ML0(1) = G0M(I,J,1)
      S0ML0(1) = S0M(I,J,1)
      UL(1,1) = UO1(IM1,J)
      UL(1,2) = UO1(I  ,J)
      UL(1,3) = VO1(I,J-1)
      UL(1,4) = VO1(I,J  )

      ULD(1,1) = VOD1(IM1,J)
      ULD(1,2) = VOD1(I  ,J)
      ULD(1,3) = UOD1(I,J-1)
      ULD(1,4) = UOD1(I,J  )

      S0ML(1) = S0M1(I,J)
      DO 250 L=2,LMIJ
      UL0(L,1) = UT(IM1,J,L)
      UL0(L,2) = UT(I  ,J,L)
      UL0(L,3) = VT(I,J-1,L)
      UL0(L,4) = VT(I,J  ,L)

      ULD0(L,1) = VTD(IM1,J,L)
      ULD0(L,2) = VTD(I  ,J,L)
      ULD0(L,3) = UTD(I,J-1,L)
      ULD0(L,4) = UTD(I,J  ,L)

      G0ML0(L) = G0M(I,J,L)
      S0ML0(L) = S0M(I,J,L)
      UL(L,1) = UL0(L,1)
      UL(L,2) = UL0(L,2)
      UL(L,3) = UL0(L,3)
      UL(L,4) = UL0(L,4)

      ULD(L,1) = ULD0(L,1)
      ULD(L,2) = ULD0(L,2)
      ULD(L,3) = ULD0(L,3)
      ULD(L,4) = ULD0(L,4)

      S0ML(L) = S0ML0(L)
  250 CONTINUE
      G0ML(1) = G0M1(I,J,1)
      DO L=2,MIN(LSRPD,LMIJ)
        G0ML(L) = G0M1(I,J,L)  ;  END DO
      DO L=LSRPD+1,LMIJ
        G0ML(L) = G0ML0(L)  ;  END DO
#ifdef TRACERS_OCEAN
      DO N=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
        TRML1(N) = TRMO1(N,I,J)
        DO L=1,LMIJ
          TRML(L,N) = TRMO(I,J,L,N)
C****
          if (entry%t_qlimit) then
            TRML(L,N) = Max (0d0,TRML(L,N))
          end if
        END DO
      END DO
#endif

C**** Calculate surface heat and salt flux: dE(W/m^2);dS,dTr(kg/m^2/s)
      DELTAE = (G0ML0(1)-G0ML(1))*BYDXYPO(J)*BYDTS
      DELTAS = (S0ML0(1)-S0ML(1))*BYDXYPO(J)*BYDTS
#ifdef TRACERS_OCEAN
      DELTATR(:) = (TRML(1,:)-TRML1(:))*BYDXYPO(J)*BYDTS
#endif
C****
C**** Vertical mixing dependent on KPP boundary layer scheme
C****
  500 IF(LMIJ.LE.1)  GO TO 725
C**** Z0 = OGEOZ/GRAV+HOCEAN Ocean depth (m)
C****      scale depths using Z0/ZE(LMIJ)

      ZSCALE = (OGEOZ(I,J) + GRAV*HOCEAN(I,J))/(ZE(LMIJ)*GRAV)
c       OIJ(I,J,IJ_OGEOZ)=OIJ(I,J,IJ_OGEOZ)+OGEOZ(I,J)
      zgrid(0) = epsln
      hwide(0) = epsln
      byhwide(0) = 0.
      POE(0) = 0.
      DO L=1,LMIJ
        POE(L) = POE(L-1) + GRAV*MO(I,J,L)
      ENDDO
      PO(1) = 5d-1*MO(I,J,1)*GRAV
      DO L=1,LMO-1
         zgrid(L) = -0.5*ZSCALE*(ZE(L-1) + ZE(L)) ! tracer level
         hwide(L) = zgrid(L-1) - zgrid(L)       ! tracer level
         byhwide(L) = 1d0/hwide(L)
         PO(L+1) = PO(L) + 5d-1*GRAV*(MO(I,J,L)+MO(I,J,L+1))
      END DO
      zgrid(LMO)   = -0.5*ZSCALE*(ZE(LMO-1) + ZE(LMO)) ! tracer level
      hwide(LMO)   = zgrid(LMO-1) - zgrid(LMO) ! tracer level
      byhwide(LMO) = 1d0/hwide(LMO)
      zgrid(LMO+1) = -ZE(LMO)*ZSCALE
      hwide(LMO+1) = epsln
      byhwide(LMO+1) = 0.

C**** If swfrac is not used add solar from next two levels to Bo
C      DELTAE = DELTAE + FSR(2)*DELTASR
C**** else remove solar from Bo
      DELTAE = DELTAE - (1d0-FSR(2))*DELTASR

C**** Start iteration for diffusivities
      MML = MMLT   ! start with pre-source mass
      BYMML = BYMMLT

      if(use_tdiss==1) then
        do l=1,lmij
          g(l) = g0ml(l)*bymml(l)
          s(l) = s0ml(l)*bymml(l)
          m1d(l) = mml(l)*bydxypo(j)
        enddo
        call get_kvtdiss(lmij,g,s,m1d,po,tdiss(i,j),tdiss_n(i,j),
     &       kvtdiss,
     &       .false.) ! debug ?
      endif

      ITER = 0
      HBL = 0
  510 ITER =ITER + 1
      HBLP = HBL
      IF (ITER.eq.2) THEN
        MML = MML0      ! continue with post-source mass
        BYMML = BYMML0
      END IF
C**** Initiallise fields
      Shsq=0.
      dVsq=0.
      dbloc=0.
      dbsfc=0.
      Ritop=0.

C**** velocity shears
      DO L=1,LMIJ-1
        DO K=1,KMUV
          Shsq(L) = Shsq(L) + RAVM(K)*(UL(L,K)-UL(L+1,K))**2 !interface
          dVsq(L) = dVsq(L) + RAVM(K)*(UL(1,K)-UL(L  ,K))**2 !tracer pts
         END DO
      END DO
      DO K=1,KMUV
        dVsq(LMIJ)=dVsq(LMIJ)+RAVM(K)*(UL(1,K)-UL(LMIJ,K))**2
      END DO

C****    density related quantities
C****
C****    density of surface layer                                (kg/m3)
C****          rho   = rho{G(1),S(1),PO(1)}
C****    local buoyancy gradient at km interfaces:                (m/s2)
C****          dbloc = g/rho{k+1,k+1} * [ rho{k,k+1}-rho{k+1,k+1} ]
C****    buoyancy difference with respect to "zref", the surface  (m/s2)
C****          dbsfc = g/rho{k,k} * [ rho{k,k}-rho{1,k} ] (CORRECT)
C****    thermal expansion coefficient without 1/rho factor    (kg/m3/C)
C****          talph= d(rho{k,k})/d(T(k))
C****    salt expansion coefficient without 1/rho factor     (kg/m3/PSU)
C****          sbeta = d(rho{k,k})/d(S(k))

      DO L=1,LMIJ
         G(L)     = G0ML(L)*BYMML(L)
         S(L)     = S0ML(L)*BYMML(L)
         BYRHO(L) = VOLGSP(G(L),S(L),PO(L))
         RHO(L)   = 1d0/BYRHO(L)
         IF (L.ge.2) THEN
           RHOM = 1d0/VOLGSP(G(L-1),S(L-1),PO(L)) ! L-1 density wrt P_L
           RHO1 = 1d0/VOLGSP(G(1),S(1),PO(L)) ! surf density wrt P_L
           dbsfc(L)   = GRAV*(1d0-RHO1*BYRHO(L))
           dbloc(L-1) = GRAV*(1d0-RHOM*BYRHO(L))
C**** numerator of bulk richardson number on grid levels
           Ritop(L) = (zgrid(1)-zgrid(L)) * dbsfc(L)
         END IF
      END DO
C**** find mld, the mixed layer depth 
      ptdd=0.03d0       ! kg/m^3  pot.dens.diff criterion
      mld=-zgrid(lmij)
      kmld=lmij
      do l=1,lmij
c       G(L)=G0ML(L)*BYMML(L)
c       S(L)=S0ML(L)*BYMML(L)
        ptd(l) = 1d0/VOLGS(G(L),S(L))-1000d0
        if (abs(ptd(l)-ptd(1)).gt.ptdd) then
          ptdm=ptd(1)+sign(ptdd,ptd(l)-ptd(1))
          mld=-zgrid(l-1)+(zgrid(l-1)-zgrid(l))*
     &        (ptdm-ptd(l-1))/(ptd(l)-ptd(l-1)+1.d-30)
          kmld=l ! kmld is the layer number immediately below the mld
          exit
        endif
      end do
      talpha(1) =  ALPHAGSP(G(1),S(1),PO(1)) ! <0
      sbeta(1)  =   BETAGSP(G(1),S(1),PO(1))*1d3

C**** surface wind turbulent friction speed (m/s) = sqrt( |tau_0|/rho )
      Ustar = SQRT(U2rho*BYRHO(1))

C**** surface buoyancy forcing turbulent and solar (m^2/s^3)
C****    Bo    = -g*(talpha * wtsfc - sbeta * wssfc)/rho
C****    Bosol = -g*(talpha * wtsol                )/rho
C**** Bo includes all buoyancy and heat forcing
      BYSHC = 1d0/SHCGS(G(1),S(1))
      Bo    = - GRAV*BYRHO(1)**2 *(
     *       sbeta(1)*DELTAS + talpha(1)*BYSHC*DELTAE -
     *     ( sbeta(1)*S(1)   + talpha(1)*BYSHC*G(1) )*DELTAM)
      Bosol = - GRAV*BYRHO(1)**2 * talpha(1)*BYSHC*DELTASR

C**** Double diffusive option (ddmix) needs alphaDT,betaDS
C**** alphaDT  = mean -talpha * delta(temp.) at interfaces (kg/m3)
C**** betaDS   = mean sbeta  * delta(salt)  at interfaces (kg/m3)

      if (LDD) then
        TO(1)     =    TEMGSP(G(1),S(1),PO(1))
        do L=2,LMIJ
          TO(L)     =   TEMGSP(G(L),S(L),PO(L))
          talpha(L)= ALPHAGSP(G(L),S(L),PO(L)) ! <0
          sbeta(L) =  BETAGSP(G(L),S(L),PO(L))*1d3
          alphaDT(L-1) = - 5d-1 * (talpha(L-1) + talpha(L))
     $         * (TO(L-1) - TO(L))
          betaDS (L-1) = - 5d-1 * (sbeta (L-1) + sbeta(L))
     $         * (S(L-1) - S(L))
        end do

        do l=1,lmij-1
          alphaDG(l) = rho(l+1) - 1d0/VOLGSP(G(L),S(L+1),PO(L+1))
          if(abs(g(l+1)-g(l))/4185d0.gt.1d-6) then
            galpha(l) = alphaDG(l)/(g(l+1)-g(l))
            if(galpha(l).eq.0d0) then
              write(6,*) 'galpha == 0',i,j,l,g(l+1)-g(l),
     &             g(l:l+1)/4185d0,s(l:l+1)
            endif
          else
            galpha(l) = -1.
          endif
          betaDS(l) = rho(l+1) - 1d0/VOLGSP(G(L+1),S(L),PO(L+1))
          if(abs(s(l+1)-s(l)).gt.1d-9) then
            sbeta(l) = betaDS(l)/(s(l+1)-s(l))
          else
            sbeta(l) = 1000.
          endif
          rhoe = .5d0*(rho(l)+rho(l+1))
          galpha(l) = galpha(l)/rhoe
          sbeta(l) = sbeta(l)/rhoe
        enddo

      end if

#ifndef OCN_GISS_TURB
      CALL KPPMIX(LDD,ZE,zgrid,hwide,LMIJ,Shsq,dVsq,Ustar,Bo
     *     ,Bosol ,alphaDT,betaDS,dbloc,Ritop,Coriol,byhwide,
     *     AKVM,AKVS,AKVG,GHAT,HBL,KBL)
      akvc=akvs
#else
      call bldepth(ZE,zgrid,byhwide,LMIJ,dVsq,Ustar,Bo,Bosol
     *     ,dbloc,Ritop
     *     ,rib,HBL,KBL,bf)
c     gawt0=GRAV*BYRHO(1)**2*talpha(1)*BYSHC*DELTAE
c     gbws0=-GRAV*BYRHO(1)**2*sbeta(1)*(DELTAS-S0ML0(1)*BYMML(1)*DELTAM)
cc    gawt0=GRAV*BYRHO(1)**2*talpha(1)*BYSHC*(DELTAE-G(1)*DELTAM)
cc    gbws0=GRAV*BYRHO(1)**2*sbeta(1)*(-DELTAS+S(1)*DELTAM)
c-c     omfrac=(bf-Bo)/(Bosol+1d-20)
c-c     wt0=-BYRHO(1)*(BYSHC*(DELTAE+omfrac*DELTASR)-
c-c    *              (BYSHC*G(1))*DELTAM)
c-c     ws0=-BYRHO(1)*(DELTAS-S(1)*DELTAM)
C if j<jm: ul(:,1)=u(i-1,j) ul(:,2)=u(i,j) ul(:,3)=v(i,j-1) ul(:,4)=v(i,j)
C          uld(:,1)=vd(i-1,j) uld(:,2)=vd(i,j) uld(:,3)=ud(i,j-1) uld(:,4)=ud(i,j)
C ud is u interpolated to the v point and vd is v interpolated to the u point.

      strait=0   ! not in strait area (for most of the oceans)
      l=lmij
      ! at tracer grid:
      u2b=0.d0
      do k=1,kmuv
        u2b=u2b+ravm(k)*ul(l,k)**2
      end do
      ut2=ut2a(i,j)
      ustarb2=cd*(ut2*(u2b+ut2))**.5d0
      exy=exya(i,j)
      rhobot(i,j)=rho(l)
      ! at velocity grid:
      if(j.lt.jm) then
        ip1=i+1
        if(i.eq.im) ip1=1
        uob=ul(l,2)
        vob=ul(l,4)
        uodb=uld(l,4)
        vodb=uld(l,2)
        u2bx=uob**2+vodb**2
        u2by=vob**2+uodb**2
        ut2x=0.5D0*(ut2+ut2a(ip1,j))
        ut2y=0.5D0*(ut2+ut2a(i,j+1))
        taubx(i,j)=cd*(u2bx+ut2x)**.5D0*uob
        tauby(i,j)=cd*(u2by+ut2y)**.5D0*vob
      endif

      do l=1,lmij-1
         e(l)=otke(l,i,j)
      end do

      call gissmix(
      ! in:
     &    lmij,ze,zgrid,dbloc,Shsq,alphaDT,betaDS,alphaDG,rho
     &   ,ustarb2,exy,Coriol,hbl,talpha,sbeta,galpha,strait,kbl,bf,ustar
     &   ,i,j
      ! out:
     &   ,ri,rrho,bv2,akvm,akvg,akvs,akvc,e)

      buoy=0.
      do l=1,lmij-1
         otke(l,i,j)=e(l)
         buoy(l)=-(akvg(l)-rrho(l)*akvs(l))/(1-rrho(l))*bv2(l)
      end do
#endif
#ifdef OCN_GISS_SM
C**** velocity at the cell center
      uc=0.d0; vc=0.d0;
      do l=1,lmij
        if(j.lt.jm) then
          uc(l)=.5d0*(ul(l,1)+ul(l,2))
          vc(l)=.5d0*(ul(l,3)+ul(l,4))
        else
          uc(l)=ul(l,im+1)
          vc(l)=ul(l,im+2)
        endif
      end do
      if(j.lt.jm) then
        rxl(:)=rx(i,j,:)
        ryl(:)=ry(i,j,:)
        gxl(:)=gx(i,j,:)
        gyl(:)=gy(i,j,:)
        sxl(:)=sx(i,j,:)
        syl(:)=sy(i,j,:)
      else
        rxl(:)=rx(im,jm,:)
        ryl(:)=rx(ivnp,jm,:)
        gxl(:)=gx(im,jm,:)
        gyl(:)=gx(ivnp,jm,:)
        sxl(:)=sx(im,jm,:)
        syl(:)=sx(ivnp,jm,:)
      endif

C-- time-averaging of u,v
      au_sm(i,j,:)=wta*au_sm(i,j,:)+(1.d0-wta)*uc(:)
      av_sm(i,j,:)=wta*av_sm(i,j,:)+(1.d0-wta)*vc(:)
      rx_sm(i,j,:)=wta*rx_sm(i,j,:)+(1.d0-wta)*rxl(:)
      ry_sm(i,j,:)=wta*ry_sm(i,j,:)+(1.d0-wta)*ryl(:)
      gx_sm(i,j,:)=wta*gx_sm(i,j,:)+(1.d0-wta)*gxl(:)
      gy_sm(i,j,:)=wta*gy_sm(i,j,:)+(1.d0-wta)*gyl(:)
      sx_sm(i,j,:)=wta*sx_sm(i,j,:)+(1.d0-wta)*sxl(:)
      sy_sm(i,j,:)=wta*sy_sm(i,j,:)+(1.d0-wta)*syl(:)
      uc(:)=au_sm(i,j,:)
      vc(:)=av_sm(i,j,:)
      rxl(:)=rx_sm(i,j,:)
      ryl(:)=ry_sm(i,j,:)
      gxl(:)=gx_sm(i,j,:)
      gyl(:)=gy_sm(i,j,:)
      sxl(:)=sx_sm(i,j,:)
      syl(:)=sy_sm(i,j,:)

      lmix=kmld
c     lmix=min(kmld,lmij-1)

      call giss_sm_mix(
      ! in:
     &    ze,zgrid,dbloc,uc,vc,rxl,ryl,gxl,gyl,sxl,syl
     &   ,grav,Coriol,lmix,lmij
      ! out:
     &   ,fvg,fvs,fvb,p4uv)

#endif

C**** Calculate non-local transport terms for each scalar
C****        ghat[sg] = kv * ghat * <w[sg]0>   (J,kg)
C**** Correct units for diffusivities (m^2/s) => (kg^2/m^4 s)
C****                            ghat (s/m^2) => (s m^4/kg^2)
      DO L=1,LMIJ-1
#ifdef ENHANCED_DEEP_MIXING
         akvg(l) = akvg(l) + kvextra(l)
         akvs(l) = akvs(l) + kvextra(l)
         akvc(l) = akvc(l) + kvextra(l)
#endif
         klen(i,j,l) = akvs(l)+kvtdiss(l)
         R = 5d-1*(RHO(L)+RHO(L+1))
         R2 = R**2
c        GHATM(L) = 0.          ! no non-local momentum transport
         DO K=1,KMUV
           GHATM(L,K) = 0. ! no non-local momentum transport
         END DO
#ifndef OCN_GISS_TURB
C**** GHAT terms must be zero for consistency with OSOURC
! why is AKV[GS]*GHAT  IF(AKVG(L)*GHAT(L) .GT. 1D0) GHAT(L)=1D0/AKVG(L)
! sometimes > 1?       IF(AKVS(L)*GHAT(L) .GT. 1D0) GHAT(L)=1D0/AKVS(L)
         GHATG(L) = AKVG(L)*GHAT(L)*DELTAE*DXYPO(J)
         GHATS(L) = AKVS(L)*GHAT(L)*(DELTAS-S0ML0(1)*BYMML(1)*DELTAM)
     *        *DXYPO(J)
#else
         GHATG(L)=0.d0
         GHATS(L)=0.d0
#endif
#ifdef TRACERS_OCEAN
#ifndef OCN_GISS_TURB
         GHATT(L,:)= AKVS(L)*GHAT(L)*(DELTATR(:)-
     *        TRML(1,:)*DELTAM*BYMML(1))*DXYPO(J)
#else
         GHATT(L,:)=0.
c        GHATT(L,:)=GHATS(L)*(DELTATR(:)-TRML(1,:)*DELTAM*BYMML(1))
c    &                      /(DELTAS-S0ML0(1)*BYMML(1)*DELTAM+1d-30)
#endif
#endif
         AKVM(L) = (AKVM(L)+KVTDISS(L))*R2
         AKVG(L) = (AKVG(L)+KVTDISS(L))*R2
         AKVS(L) = (AKVS(L)+KVTDISS(L))*R2
         AKVC(L) = (AKVC(L)+KVTDISS(L))*R2
         AKVG3D(L,I,J) = AKVG(L)
         AKVS3D(L,I,J) = AKVS(L)
         AKVC3D(L,I,J) = AKVC(L)
      END DO
      AKVG3D(0,I,J) = AKVG3D(1,I,J)
      AKVS3D(0,I,J) = AKVS3D(1,I,J)
      AKVC3D(0,I,J) = AKVC3D(1,I,J)
      AKVG3D(LMIJ,I,J) = AKVG3D(LMIJ-1,I,J)
      AKVS3D(LMIJ,I,J) = AKVS3D(LMIJ-1,I,J)
      AKVC3D(LMIJ,I,J) = AKVC3D(LMIJ-1,I,J)

      DTP4G=0.; DTP4S=0.; DTP4UV=0.
#ifdef TRACERS_OCEAN
      DTP4TR=0.
#endif
#ifdef OCN_GISS_SM
      ! at the layer middle:
      ! if j=jm, then i stays at 1
      DO L=1,LMIJ-1
         bydz=DTS/(ze(l-1)-ze(l))*MML(L)
         DTP4G(l)=(fvg(l-1)-fvg(l))*bydz
         DTP4S(l)=(fvs(l-1)-fvs(l))*bydz
c        DTP4G(l)=-DTS*dfvgdz(l)*MML(L)
c        DTP4S(l)=-DTS*dfvsdz(l)*MML(L)
c        DTP4G(l)=0.d0
c        DTP4S(l)=0.d0
c        DO K=1,KMUV
c          DTP4UV(L,K) = 0.d0
c        END DO
         DTP4G3D(L,I,J) = DTP4G(L)
         DTP4S3D(L,I,J) = DTP4S(L)
      END DO
      DTP4G3D(0,I,J) = DTP4G3D(1,I,J)
      DTP4S3D(0,I,J) = DTP4S3D(1,I,J)
      DTP4G3D(LMIJ,I,J) = DTP4G3D(LMIJ-1,I,J)
      DTP4S3D(LMIJ,I,J) = DTP4S3D(LMIJ-1,I,J)
      DO L=1,LMIJ
         DO K=1,KMUV
           if(j.lt.jm.and.k.le.2) then
             DTP4UV(l,k)=DTS*p4uv(l,1)
           else
             DTP4UV(l,k)=DTS*p4uv(l,2)
           endif
         END DO
      END DO
#endif

C**** For each field (U,G,S + TRACERS) call OVDIFF
C**** Momentum
      DO K=1,KMUV
        IF(LMUV(K).GT.1) THEN
c         CALL OVDIFF(UL(1,K),AKVM(1),GHATM,DTBYDZ,BYDZ2
c    *         ,LMUV(K),UL0(1,K))
          CALL OVDIFF(UL(1,K),AKVM(1),GHATM(1,K)
     *        ,DTP4UV(1,K),ZE(1),ZGRID(1)
     *        ,DTBYDZ,BYDZ2,LMUV(K),UL0(1,K))
        ENDIF
      END DO
C**** Enthalpy
      Call OVDIFFS (G0ML(1),AKVG(1),GHATG,DTP4G,DTBYDZ,BYDZ2,DTS,LMIJ,
     *             G0ML0(1),FLG)
C**** Salinity
      Call OVDIFFS (S0ML(1),AKVS(1),GHATS,DTP4S,DTBYDZ,BYDZ2,DTS,LMIJ,
     *              S0ML0(1),FLS)
      IF ((ITER.eq.1  .or. ABS(HBLP-HBL).gt.(ZE(KBL)-ZE(KBL-1))*0.25)
     *     .and. ITER.lt.4) GO TO 510
C**** D-grid velocities
      If (.not.QPOLE)  Then
        DO K=1,KMUV
          IF(LMUV(K).GT.1) THEN
            CALL OVDIFF(ULD(1,K),AKVM(1),GHATM(1,K)
     *       ,DTP4UV(1,K),ZE(1),ZGRID(1)
     *       ,DTBYDZ,BYDZ2,LMUV(K),ULD0(1,K))
          ENDIF
        END DO

      EndIf
#ifdef TRACERS_OCEAN
C**** Tracers are diffused after iteration and follow salinity
      if(ntrtrans.gt.1) then
        do l=1,lmij
          dtbydz_tr(l) = dts/motr(i,j,l)
        enddo
        do l=1,lmij-1
          bydz2_tr(l)  = 2d0/(motr(i,j,l)+motr(i,j,l+1))
        enddo
      else
        do l=1,lmij
          dtbydz_tr(l) = dtbydz(l)
          bydz2_tr(l) = bydz2(l)
        enddo
      endif
      GHATDUM(:) = 0.
      DTP4S(:)   = 0.  ! ????
      DO N=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
        if(entry%t_qlimit) then
          ! Modify GHATT to prevent negative tracer.  Method: apply
          ! a single multiplicative factor to the GHATT profile.
          ! Might switch to alternative method (local flux adjustments).
          minrat = 1.
          dflux = ghatt(1,n)*dts
          if(dflux.gt.0. .and. trml(1,n).lt.dflux) then
            minrat = trml(1,n)/dflux
          endif
          do l=2,lmij-1
            dflux = (ghatt(l,n)-ghatt(l-1,n))*dts
            if(dflux.gt.0. .and. trml(l,n).lt.dflux) then
              minrat = min(minrat, trml(l,n)/dflux)
            endif
          enddo
          minrat = max(minrat, 0.d0) ! in case min tracer was already < 0
          if(minrat .lt. 1.d0) then
            minrat = minrat*0.95d0 ! leave min tracer slightly > 0
            ghatt(1:lmij-1,n) = ghatt(1:lmij-1,n)*minrat
          endif
        endif
        ! diffuse
        Call OVDIFFS (TRML(1,N),AKVC(1),GHATT(1,N),DTP4TR(1,N)
     *       ,DTBYDZ_tr,BYDZ2_tr,DTS,LMIJ,TRML(1,N),FLT(1,N))
      END DO
#endif

C**** Save fluxes for adjustment of vertical gradients
      DM = DELTAM*DTBYDZ(1)
      FLG3D(1:LMIJ-1,I,J) = FLG(1:LMIJ-1)
      FLG3D(LMIJ,I,J) = 0.
      FLG3D(0,I,J) = DTS*DELTAE*DXYPO(J)
      FLS3D(1:LMIJ-1,I,J) = FLS(1:LMIJ-1)
      FLS3D(LMIJ,I,J) = 0.
      FLS3D(0,I,J) = -DTS*DELTAS*DXYPO(J)*(1-DM) +DM*S0M1(I,J)
#ifdef TRACERS_OCEAN
      FLT3D(1:LMIJ-1,:,I,J) = FLT(1:LMIJ-1,:)
      FLT3D(LMIJ,:,I,J) = 0.
      FLT3D(0,:,I,J) = -DTS*DELTATR(:)*DXYPO(J)*(1-DM) +DM*TRMO1(:,I,J)
#endif

      DZ3D(1:LMIJ,I,J) = MO(I,J,1:LMIJ)*BYRHO(1:LMIJ)

C**** Diagnostics for non-local transport and vertical diffusion
       DO L=1,LMIJ-1
         OIJL(I,J,L,IJL_KVM) = OIJL(I,J,L,IJL_KVM) + AKVM(L)
         OIJL(I,J,L,IJL_KVG) = OIJL(I,J,L,IJL_KVG) + AKVG(L)
         if(use_tdiss==1) then
           OIJL(I,J,L,IJL_KVX) = OIJL(I,J,L,IJL_KVX) + KVTDISS(L)
         endif
#ifdef OCN_GISS_TURB
         OIJL(I,J,L,IJL_KVS) = OIJL(I,J,L,IJL_KVS) + AKVS(L)
         OIJL(I,J,L,IJL_KVC) = OIJL(I,J,L,IJL_KVC) + AKVC(L)
         OIJL(I,J,L,ijl_ri) = OIJL(I,J,L,ijl_ri) + ri(L) ! Richardson number
         OIJL(I,J,L,ijl_rrho)= OIJL(I,J,L,ijl_rrho) + rrho(L) ! density ratio
         OIJL(I,J,L,ijl_bv2)= OIJL(I,J,L,ijl_bv2) + bv2(L) ! B.Vaisala freq sq
         OIJL(I,J,L,ijl_buoy)= OIJL(I,J,L,ijl_buoy) + buoy(L) ! buoyancy flux
         OIJL(I,J,L,ijl_otke)= OIJL(I,J,L,ijl_otke) + e(L) ! turbulent k.e.
#endif
#ifdef OCN_GISS_SM
         OIJL(I,J,L,ijl_fvb)= OIJL(I,J,L,ijl_fvb) + fvb(L)
#endif
         OIJL(I,J,L,IJL_WGFL)= OIJL(I,J,L,IJL_WGFL) + FLG(L) ! heat flux
         OIJL(I,J,L,IJL_WSFL)= OIJL(I,J,L,IJL_WSFL) + FLS(L) ! salt flux
#ifdef TRACERS_OCEAN
C**** vertical diffusive tracer flux
         TOIJL(I,J,L,TOIJL_WTFL,:)=TOIJL(I,J,L,TOIJL_WTFL,:)+FLT(L,:)
#endif
       END DO
C**** Also set vertical diagnostics
       DO L=1,LMIJ
CCC         OL(L,L_RHO) = OL(L,L_RHO) + RHO(L)
CCC         OL(L,L_TEMP)= OL(L,L_TEMP)+ TEMGS(G(L),S(L))
CCC         OL(L,L_SALT)= OL(L,L_SALT)+ S(L)
         OLJ(1,L,J)= OLJ(1,L,J) + RHO(L) ! L_RHO
         OLJ(2,L,J)= OLJ(2,L,J) + TEMGS(G(L),S(L)) ! L_TEMP
         OLJ(3,L,J)= OLJ(3,L,J) + S(L) ! L_SALT
       END DO
C**** Set diagnostics
       OIJ(I,J,IJ_HBL) = OIJ(I,J,IJ_HBL) + HBL ! boundary layer depth
       OIJ(I,J,ij_mld) = OIJ(I,J,ij_mld) + mld ! mixed layer depth
       OIJ(I,J,IJ_BO) = OIJ(I,J,IJ_BO) + Bo ! surface buoyancy forcing
       OIJ(I,J,IJ_BOSOL) = OIJ(I,J,IJ_BOSOL) + Bosol ! solar buoy frcg
       OIJ(I,J,IJ_USTAR) = OIJ(I,J,IJ_USTAR) + Ustar ! turb fric speed

       OIJmm(I,J,IJ_HBLmax) = max(OIJmm(I,J,IJ_HBLmax), HBL)
       OIJmm(I,J,ij_mldmax) = max(OIJmm(I,J,IJ_mldmax), mld)

       IF(KBL.gt.KPL(I,J)) KPL(I,J)=KBL ! save max. mixed layer depth
C****
C**** Update current prognostic variables
C****
CCC      IF (QPOLE) THEN
CCC        DO II=1,IM
CCC          VO(II,JM-1,1:LMUV(II))=VO(II,JM-1,1:LMUV(II))+RAMV(II)*
CCC     *         (UL(1:LMUV(II),II)-VT(II,JM-1,1:LMUV(II)))
CCC        END DO
CCC        UO(IM,JM,1:LMIJ)=UO(IM,JM,1:LMIJ) + RAMV(IM+1)*
CCC     *       (UL(1:LMIJ,IM+1)-UT(IM,JM,1:LMIJ))
CCC        UO(IVNP,JM,1:LMIJ)=UO(IVNP,JM,1:LMIJ) + RAMV(IM+2)*
CCC     *       (UL(1:LMIJ,IM+2)-UT(IVNP,JM,1:LMIJ))
CCC      ELSE
CCC        UO(IM1,J,1:LMUV(1))=UO(IM1,J,1:LMUV(1)) + RAMV(1)*
CCC     *       (UL(1:LMUV(1),1)-UT(IM1,J,1:LMUV(1)))
CCC        UO(I  ,J,1:LMUV(2))=UO(I  ,J,1:LMUV(2)) + RAMV(2)*
CCC     *       (UL(1:LMUV(2),2)-UT(I  ,J,1:LMUV(2)))
CCC        VO(I,J-1,1:LMUV(3))=VO(I,J-1,1:LMUV(3)) + RAMV(3)*
CCC     *       (UL(1:LMUV(3),3)-VT(I,J-1,1:LMUV(3)))
CCC        VO(I,J  ,1:LMUV(4))=VO(I,J  ,1:LMUV(4)) + RAMV(4)*
CCC     *       (UL(1:LMUV(4),4)-VT(I,J  ,1:LMUV(4)))
CCC      END IF

      IF(QPOLE)  THEN
         DO II=1,IM
           UKJM(1:LMUV(II),II) = RAMV(II)*(UL(1:LMUV(II),II)-
     *          VT(II,JM-1,1:LMUV(II)))
         END DO
         UKJM(1:LMIJ,IM+1)= RAMV(IM+1)*(UL(1:LMIJ,IM+1)-
     *        UT(IM,JM,1:LMIJ))
         UKJM(1:LMIJ,IM+2)= RAMV(IM+2)*(UL(1:LMIJ,IM+2)-
     *        UT(IVNP,JM,1:LMIJ))
      ELSE
        UKM(1:LMUV(1),1,I,J) = RAMV(1)*(UL(1:LMUV(1),1)
     *          -UT(IM1,J,1:LMUV(1)))
        UKM(1:LMUV(2),2,I,J) = RAMV(2)*(UL(1:LMUV(2),2)
     *          -UT(I  ,J,1:LMUV(2)))
        UKM(1:LMUV(3),3,I,J) = RAMV(3)*(UL(1:LMUV(3),3)
     *          -VT(I,J-1,1:LMUV(3)))
        UKM(1:LMUV(4),4,I,J) = RAMV(4)*(UL(1:LMUV(4),4)
     *          -VT(I,J  ,1:LMUV(4)))

        UKMD(1:LMUV(1),1,I,J) = RAMV(1)*(ULD(1:LMUV(1),1)
     *          -VTD(IM1,J,1:LMUV(1)))
        UKMD(1:LMUV(2),2,I,J) = RAMV(2)*(ULD(1:LMUV(2),2)
     *          -VTD(I  ,J,1:LMUV(2)))
        UKMD(1:LMUV(3),3,I,J) = RAMV(3)*(ULD(1:LMUV(3),3)
     *          -UTD(I,J-1,1:LMUV(3)))
        UKMD(1:LMUV(4),4,I,J) = RAMV(4)*(ULD(1:LMUV(4),4)
     *          -UTD(I,J  ,1:LMUV(4)))
      END IF

  725 IF(QPOLE)  GO TO 750
      DO L=1,LMIJ
        G0M(I,J,L) = G0ML(L)
        S0M(I,J,L) = S0ML(L)
      END DO
#ifdef TRACERS_OCEAN
      DO N = 1,tracerlist%getsize()
      DO L = 1,LMIJ
        TRMO(I,J,L,N) = TRML(L,N)
      ENDDO
      ENDDO
#endif
C**** End of I loop
  730 IM1=I
      I=I+1
      IF(I.LE.IM)  GO TO 210
      GO TO 790
C****
C**** Load the prognostic variables of potential enthalpy and
C**** salt from the column arrays at the poles
C****
  750 DO L=1,LMIJ
      G0M(1,JM,L)  = G0ML(L)
      S0M(1,JM,L)  = S0ML(L)
#ifdef TRACERS_OCEAN
      TRMO(1,JM,L,:) = TRML(L,:)
#endif
      END DO
C**** End of outside J loop
  790 CONTINUE

c#ifdef OCN_GISS_TURB
c     if(have_north_pole) then
c     ! for the north pole, in subroutine obdrag only (not active) 
c       unp=0.d0
c       vnp=0.d0
c       do i=1,im
c         unp=unp-sinic(i)*tauby(i,jm-1)
c         vnp=vnp+cosic(i)*tauby(i,jm-1)
c       end do
c       taubx(im,jm)  =unp*2/im
c       taubx(ivnp,jm)=vnp*2/im
c     endif
c#endif

C**** Update velocities outside parallel region
      call halo_update_block (grid, UKM, from=north)
      call halo_update_block (grid, UKMD, from=north)

C**** North pole
c      if (have_north_pole) then
c      DO I=1,IM
c        VO(I,JM-1,1:LMV(I,JM-1))=VO(I,JM-1,1:LMV(I,JM-1))+
c     *       UKJM(1:LMV(I,JM-1),I)
c      END DO
c      UO(IM,JM,1:LMU(1,JM))=UO(IM,JM,1:LMU(1,JM))+UKJM(1:LMU(1,JM),IM+1)
c      UO(IVNP,JM,1:LMU(1,JM)) = UO(IVNP,JM,1:LMU(1,JM)) +
c     +                          UKJM(1:LMU(1,JM),IM+2)
c      end if

C**** Everywhere else
      DO J=j_0s,j_1s
        IM1=IM
        DO I=1,IM
          UO(IM1,J,1:LMU(IM1,J)) = UO(IM1,J,1:LMU(IM1,J)) +
     +                             UKM(1:LMU(IM1,J),1,I,J)
          UO(I  ,J,1:LMU(I  ,J)) = UO(I  ,J,1:LMU(I  ,J)) +
     +                             UKM(1:LMU(I  ,J),2,I,J)
          VO(I,J-1,1:LMV(I,J-1)) = VO(I,J-1,1:LMV(I,J-1)) +
     +                             UKM(1:LMV(I,J-1),3,I,J)
          VO(I,J  ,1:LMV(I,J  )) = VO(I,J  ,1:LMV(I,J  )) +
     +                             UKM(1:LMV(I,J  ),4,I,J)
#ifndef DISABLE_KPP_DGRID_MIXING
          VOD(IM1,J,1:LMU(IM1,J)) = VOD(IM1,J,1:LMU(IM1,J)) +
     +                             UKMD(1:LMU(IM1,J),1,I,J)
          VOD(I  ,J,1:LMU(I  ,J)) = VOD(I  ,J,1:LMU(I  ,J)) +
     +                             UKMD(1:LMU(I  ,J),2,I,J)
          UOD(I,J-1,1:LMV(I,J-1)) = UOD(I,J-1,1:LMV(I,J-1)) +
     +                             UKMD(1:LMV(I,J-1),3,I,J)
          UOD(I,J  ,1:LMV(I,J  )) = UOD(I,J  ,1:LMV(I,J  )) +
     +                             UKMD(1:LMV(I,J  ),4,I,J)
#endif
          IM1=I
        END DO
      END DO
      if (.not.have_north_pole) then
        j=j_1s+1
        DO I=1,IM
          VO(I,J-1,1:LMV(I,J-1)) = VO(I,J-1,1:LMV(I,J-1)) +
     *                             UKM(1:LMV(I,J-1),3,I,J)
#ifndef DISABLE_KPP_DGRID_MIXING
          UOD(I,J-1,1:LMV(I,J-1)) = UOD(I,J-1,1:LMV(I,J-1)) +
     *                             UKMD(1:LMV(I,J-1),3,I,J)
#endif
        END DO
      end if

C**** sum global mean diagnostics
c     DO J=j_0s,j_1
c       DO L=1,LMO
c         OL(L,L_RHO) = OL(L,L_RHO) + OLJ(1,L,J)
c         OL(L,L_TEMP)= OL(L,L_TEMP)+ OLJ(2,L,J)
c         OL(L,L_SALT)= OL(L,L_SALT)+ OLJ(3,L,J)
c       END DO
c     END DO

      call globalsum(grid,OLJ(1,:,:),OLtemp(:,L_RHO)  ,jband=(/1,jm/) )
      call globalsum(grid,OLJ(2,:,:),OLtemp(:,L_TEMP) ,jband=(/1,jm/) )
      call globalsum(grid,OLJ(3,:,:),OLtemp(:,L_SALT) ,jband=(/1,jm/) )
      OL(:,L_RHO ) = OL(:,L_RHO ) + OLtemp(:,L_RHO)
      OL(:,L_TEMP) = OL(:,L_TEMP) + OLtemp(:,L_TEMP)
      OL(:,L_SALT) = OL(:,L_SALT) + OLtemp(:,L_SALT)


      if(adjust_zslope_using_flux) then
C**** Vertical Gradients of scalars
C**** Implicitly apply interpolated KV to linear profile
C**** Surface tracer + mass fluxes included (no solar flux)
C**** Note that FL[GS] are upward fluxes.
#ifdef TRACERS_OCEAN
        if(ntrtrans.gt.1) then
          call stop_model('ocnkpp: not yet implemented',255)
        endif
#endif
        do j=j_0,j_1
        do ib=1,nbyzm(j,1)
        do i=i1yzm(ib,j,1),i2yzm(ib,j,1)
          lmij = lmm(i,j)
          dtbydz(1:lmij) = dts/mo(i,j,1:lmij)
          DO L=1,LMIJ
            DTBYDZ2 = 6d0*DTBYDZ(L)**2*BYDTS
            GZMO(I,J,L)=(GZMO(I,J,L)+3d0*(FLG3D(L-1,I,J)+FLG3D(L,I,J)))
     *           /(1d0+DTBYDZ2*(AKVG3D(L-1,I,J)+AKVG3D(L,I,J)))
            SZMO(I,J,L)=(SZMO(I,J,L)+3d0*(FLS3D(L-1,I,J)+FLS3D(L,I,J)))
     *           /(1d0+DTBYDZ2*(AKVS3D(L-1,I,J)+AKVS3D(L,I,J)))
#ifdef TRACERS_OCEAN
            TZMO(I,J,L,:)=(TZMO(I,J,L,:)
     &           +3d0*(FLT3D(L-1,:,I,J)+FLT3D(L,:,I,J)))
     &           /(1d0+DTBYDZ2*(AKVC3D(L-1,I,J)+AKVC3D(L,I,J)))
#endif
          END DO
        enddo
        enddo
        enddo
      elseif(relax_subgrid_zprofile) then
        do j=j_0,j_1
        do n=1,nbyzm(j,1)
        do i=i1yzm(n,j,1),i2yzm(n,j,1)
          l = lmm(i,j)
          klen(i,j,l) = klen(i,j,l-1)
        enddo
        enddo
        enddo
        do l=lmo-1,2,-1
        do j=j_0,j_1
        do n=1,nbyzm(j,l+1)
        do i=i1yzm(n,j,l+1),i2yzm(n,j,l+1)
          klen(i,j,l) = .5d0*(klen(i,j,l)+klen(i,j,l-1))
        enddo
        enddo
        enddo
        enddo
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          MA(I,J,L) = MO(I,J,L)*DXYPO(J)
          !klen(i,j,l) = exp(-12.*klen(i,j,l)*dts/(dzo(l)*dzo(l)))
          klen(i,j,l) = 1./(1.+12.*klen(i,j,l)*dts/(dzo(l)*dzo(l)))
        enddo
        enddo
        enddo
        enddo
        call relax_zmoms(ma,g0m,klen,gzmo,gzzmo)
        call relax_zmoms(ma,s0m,klen,szmo,szzmo)
#ifdef TRACERS_OCEAN
        if(ntrtrans.gt.1) then
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          ma(i,j,l) = motr(i,j,l)*dxypo(j)
        enddo
        enddo
        enddo
        enddo
        endif
        do n=1,tracerlist%getsize()
          call relax_zmoms(ma,trmo(1,j_0h,1,n),
     &         klen,tzmo(1,j_0h,1,n),tzzmo(1,j_0h,1,n))
        enddo
#endif
      elseif(mix_tripled_resolution) then
#ifdef TRACERS_OCEAN
        if(ntrtrans.gt.1) then
          call stop_model('ocnkpp: not yet implemented',255)
        endif
#endif
        do j=j_0,j_1
        do ib=1,nbyzm(j,1)
        do i=i1yzm(ib,j,1),i2yzm(ib,j,1)
          lmij = lmm(i,j)
          !akvg(1:lmij) = akvg3d(1:lmij,i,j)
          akvg(1:lmij) = klen(i,j,1:lmij) ! m2/s units!!!!
          mml(1:lmij) = mo(i,j,1:lmij)*dxypo(j)
c
          g0ml(1:lmij) = gsave3d(i,j,1:lmij)/mml(1:lmij)
          gzml(1:lmij) = gzmo(i,j,1:lmij)/mml(1:lmij)
          gzzml(1:lmij) = gzzmo(i,j,1:lmij)/mml(1:lmij)
          call diffuse_moms(dz3d(1,i,j),akvg,g0ml,gzml,gzzml,lmij,dts)
          gzmo(i,j,1:lmij) = gzml(1:lmij)*mml(1:lmij)
          gzzmo(i,j,1:lmij) = gzzml(1:lmij)*mml(1:lmij)
c
          s0ml(1:lmij) = ssave3d(i,j,1:lmij)/mml(1:lmij)
          szml(1:lmij) = szmo(i,j,1:lmij)/mml(1:lmij)
          szzml(1:lmij) = szzmo(i,j,1:lmij)/mml(1:lmij)
          call diffuse_moms(dz3d(1,i,j),akvg,s0ml,szml,szzml,lmij,dts)
          szmo(i,j,1:lmij) = szml(1:lmij)*mml(1:lmij)
          szzmo(i,j,1:lmij) = szzml(1:lmij)*mml(1:lmij)
#ifdef TRACERS_OCEAN
          do n=1,tracerlist%getsize()
            trml(1:lmij,1) = trsave3d(i,j,1:lmij,n)/mml(1:lmij)
            tzml(1:lmij,1) = tzmo(i,j,1:lmij,n)/mml(1:lmij)
            tzzml(1:lmij,1) = tzzmo(i,j,1:lmij,n)/mml(1:lmij)
            call diffuse_moms(dz3d(1,i,j),akvg,trml,tzml,tzzml,lmij,dts)
            tzmo(i,j,1:lmij,n) = tzml(1:lmij,1)*mml(1:lmij)
            tzzmo(i,j,1:lmij,n) = tzzml(1:lmij,1)*mml(1:lmij)
          enddo
#endif
        enddo
        enddo
        enddo
      endif


C****
C**** Vertically mix horizontal gradients as though they were scalars
C****
      ghatdum = 0.
      do j=j_0s,j_1s ! gradients always set to zero in polar boxes
      do ib=1,nbyzm(j,1)
      do i=i1yzm(ib,j,1),i2yzm(ib,j,1)
        lmij = lmm(i,j)
        dtbydz(1:lmij) = dts/mo(i,j,1:lmij)
        bydz2(1:lmij-1)= 2d0/(mo(i,j,2:lmij)+mo(i,j,1:lmij-1))
        akvg(1:lmij-1) = akvg3d(1:lmij-1,i,j)
        akvg(lmij) = 0.
        akvs(1:lmij-1) = akvs3d(1:lmij-1,i,j)
        akvs(lmij) = 0.
#ifdef OCN_GISS_SM
        DTP4G(1:lmij-1) = DTP4G3D(1:lmij-1,i,j)
        DTP4G(lmij) = 0.
        DTP4S(1:lmij-1) = DTP4S3D(1:lmij-1,i,j)
        DTP4S(lmij) = 0.
#endif
        gxml(1:lmij) = gxmo(i,j,1:lmij)
        gyml(1:lmij) = gymo(i,j,1:lmij)
        sxml(1:lmij) = sxmo(i,j,1:lmij)
        syml(1:lmij) = symo(i,j,1:lmij)
        Call OVDIFFS (GXML(1),AKVG(1),GHATDUM,DTP4G,DTBYDZ,BYDZ2
     *                ,DTS,LMIJ,GXML(1),FLDUM)
        Call OVDIFFS (GYML(1),AKVG(1),GHATDUM,DTP4G,DTBYDZ,BYDZ2
     *                ,DTS,LMIJ,GYML(1),FLDUM)
        Call OVDIFFS (SXML(1),AKVS(1),GHATDUM,DTP4S,DTBYDZ,BYDZ2
     *                ,DTS,LMIJ,SXML(1),FLDUM)
        Call OVDIFFS (SYML(1),AKVS(1),GHATDUM,DTP4S,DTBYDZ,BYDZ2
     *                ,DTS,LMIJ,SYML(1),FLDUM)
        gxmo(i,j,1:lmij) = gxml(1:lmij)
        gymo(i,j,1:lmij) = gyml(1:lmij)
        sxmo(i,j,1:lmij) = sxml(1:lmij)
        symo(i,j,1:lmij) = syml(1:lmij)
        if(use_qus==1) then
        gxxml(1:lmij) = gxxmo(i,j,1:lmij)
        gyyml(1:lmij) = gyymo(i,j,1:lmij)
        gxyml(1:lmij) = gxymo(i,j,1:lmij)
        sxxml(1:lmij) = sxxmo(i,j,1:lmij)
        syyml(1:lmij) = syymo(i,j,1:lmij)
        sxyml(1:lmij) = sxymo(i,j,1:lmij)
        Call OVDIFFS (GXXML(1),AKVG(1),GHATDUM,DTP4G,DTBYDZ,BYDZ2
     *                ,DTS,LMIJ,GXXML(1),FLDUM)
        Call OVDIFFS (GYYML(1),AKVG(1),GHATDUM,DTP4G,DTBYDZ,BYDZ2
     *                ,DTS,LMIJ,GYYML(1),FLDUM)
        Call OVDIFFS (GXYML(1),AKVG(1),GHATDUM,DTP4G,DTBYDZ,BYDZ2
     *                ,DTS,LMIJ,GXYML(1),FLDUM)
        Call OVDIFFS (SXXML(1),AKVS(1),GHATDUM,DTP4S,DTBYDZ,BYDZ2
     *                ,DTS,LMIJ,SXXML(1),FLDUM)
        Call OVDIFFS (SYYML(1),AKVS(1),GHATDUM,DTP4S,DTBYDZ,BYDZ2
     *                ,DTS,LMIJ,SYYML(1),FLDUM)
        Call OVDIFFS (SXYML(1),AKVS(1),GHATDUM,DTP4S,DTBYDZ,BYDZ2
     *                ,DTS,LMIJ,SXYML(1),FLDUM)
        gxxmo(i,j,1:lmij) = gxxml(1:lmij)
        gyymo(i,j,1:lmij) = gyyml(1:lmij)
        gxymo(i,j,1:lmij) = gxyml(1:lmij)
        sxxmo(i,j,1:lmij) = sxxml(1:lmij)
        syymo(i,j,1:lmij) = syyml(1:lmij)
        sxymo(i,j,1:lmij) = sxyml(1:lmij)
        endif ! using qus
#ifdef TRACERS_OCEAN
        akvc(1:lmij-1) = akvc3d(1:lmij-1,i,j)
        akvc(lmij) = 0.
        do n=1,tracerlist%getsize()
          txml(1:lmij) = txmo(i,j,1:lmij,n)
          tyml(1:lmij) = tymo(i,j,1:lmij,n)
          Call OVDIFFS ( TXML(1),AKVC(1),GHATDUM,DTP4S,DTBYDZ,BYDZ2,
     *         DTS,LMIJ, TXML(1),FLDUM)
          Call OVDIFFS ( TYML(1),AKVC(1),GHATDUM,DTP4S,DTBYDZ,BYDZ2,
     *         DTS,LMIJ, TYML(1),FLDUM)
          txmo(i,j,1:lmij,n) = txml(1:lmij)
          tymo(i,j,1:lmij,n) = tyml(1:lmij)
          if(use_qus==1) then
          txxml(1:lmij) = txxmo(i,j,1:lmij,n)
          tyyml(1:lmij) = tyymo(i,j,1:lmij,n)
          txyml(1:lmij) = txymo(i,j,1:lmij,n)
          Call OVDIFFS (TXXML(1),AKVC(1),GHATDUM,DTP4S,DTBYDZ,BYDZ2,
     *         DTS,LMIJ,TXXML(1),FLDUM)
          Call OVDIFFS (TYYML(1),AKVC(1),GHATDUM,DTP4S,DTBYDZ,BYDZ2,
     *         DTS,LMIJ,TYYML(1),FLDUM)
          Call OVDIFFS (TXYML(1),AKVC(1),GHATDUM,DTP4S,DTBYDZ,BYDZ2,
     *         DTS,LMIJ,TXYML(1),FLDUM)
          txxmo(i,j,1:lmij,n) = txxml(1:lmij)
          tyymo(i,j,1:lmij,n) = tyyml(1:lmij)
          txymo(i,j,1:lmij,n) = txyml(1:lmij)
          endif ! using qus
        enddo
#endif
      enddo
      enddo
      enddo

      if(extra_slope_limitations) then
        do l=1,lmo
        do j=j_0s,j_1s
        do ib=1,nbyzm(j,l)
        do i=i1yzm(ib,j,l),i2yzm(ib,j,l)
C****
C**** Recreate main prognostic variables of potential enthalpy and
C**** salt from those on the quarter boxes
C****
          NSIGG = EXPONENT(G0M(I,J,L)) -2 - 42
          CALL REDUCE_FIG(NSIGG,GXMO(I,J,L))
          CALL REDUCE_FIG(NSIGG,GYMO(I,J,L))
          NSIGS = EXPONENT(S0M(I,J,L)) - 2 - 42 + 4
          CALL REDUCE_FIG(NSIGS,SXMO(I,J,L))
          CALL REDUCE_FIG(NSIGS,SYMO(I,J,L))
C**** limit salinity gradients
          TXY = abs(SXMO(I,J,L)) + abs(SYMO(I,J,L))
          if ( TXY > S0M(I,J,L) ) then
            SXMO(I,J,L)=SXMO(I,J,L)*(S0M(I,J,L)/(TXY + tiny(TXY)))
            SYMO(I,J,L)=SYMO(I,J,L)*(S0M(I,J,L)/(TXY + tiny(TXY)))
          end if
          if ( abs(SZMO(I,J,L)) > S0M(I,J,L) )
     *         SZMO(I,J,L) = sign(S0M(I,J,L),SZMO(I,J,L)+0d0)
#ifdef TRACERS_OCEAN
          DO N = 1,tracerlist%getsize()
            entry=>tracerlist%at(n)
            NSIGT = EXPONENT(TRMO(I,J,L,N)) - 2 - 42
            CALL REDUCE_FIG(NSIGT,TXMO(I,J,L,N))
            CALL REDUCE_FIG(NSIGT,TYMO(I,J,L,N))
C****
            if (entry%t_qlimit) then ! limit gradients
              TXY = abs(TXMO(I,J,L,N)) + abs(TYMO(I,J,L,N))
              if ( TXY > TRMO(I,J,L,N) ) then
                TXMO(I,J,L,N) = TXMO(I,J,L,N)
     *               *( TRMO(I,J,L,N)/(TXY + tiny(TXY)) )
                TYMO(I,J,L,N) = TYMO(I,J,L,N)
     *               *( TRMO(I,J,L,N)/(TXY + tiny(TXY)) )
              end if
              if ( abs(TZMO(I,J,L,N)) > TRMO(I,J,L,N) )
     *             TZMO(I,J,L,N) = sign(TRMO(I,J,L,N),TZMO(I,J,L,N)+0d0)
            end if
C****
          END DO
#endif
        enddo
        enddo
        enddo
        enddo
      endif ! extra_slope_limitations
C****

      RETURN
      END SUBROUTINE OCONV

      subroutine get_kvtdiss(lm,g,s,m,p,tpow,tpow_n,kv,do_debug)
!@sum get_kvtdiss Deliberately simplified calculation of tidally
!@+   generated diffusivity kv from a semi-interactive estimate
!@+   of column-integrated tidal dissipation.
!@+   This version employs a bulk (nonlocal) relationship
!@+   between dissipation and diffusivity, in contrast to
!@+   other modelE codes which set kv = dissipation/N2
!@+   locally at each depth with dissipation independent
!@+   of local N2.  Here, the rate of column potential energy
!@+   gain from the action of a trial profile of kv is
!@+   calculated; kv is then multiplicatively adjusted
!@+   to match the column-integrated dissipation rate.  The
!@+   shape of the kv profile is an exponential decay upward
!@+   from the local seafloor with a vertical scale inversely
!@+   proportional to column N.  The reference dissipation rate
!@+   in each column is taken from an offline model and is
!@+   rescaled by the ratio N/N_offline where N is the
!@+   near-seafloor B-V frequency.
!@+   Note that output kv is an effective buoyancy diffusivity,
!@+   not a diffusivity for heat, salt, or tracers individually.
!@+   Consistent with the host KPP code, this version assumes equality
!@+   of heat and salt diffusivity when evaluating the action of the
!@+   trial kv.
      use constant, only : grav
      use kpp_com, only : eff => tdiss_eff
      implicit none
!@var lm number of layers in local column
      integer :: lm
!@var tpow column-integrated tidal dissipation (W/m2) from offline model
!@var tpow_n B-V frequency (1/s) used to compute tpow in offline model
      real*8 :: tpow,tpow_n
!@var g specific potential enthalphy (J/kg)
!@var s salinity (kg/kg)
!@var m layer mass (kg/m2)
!@var p pressure (Pa)
      real*8, dimension(lm) :: g,s,m,p
!@var kv output tidally induced diffusivity (m2/s)
      real*8, dimension(lm) :: kv
!@var do_debug whether to output debugging info
      logical :: do_debug
c
! Empirical parameters and tuning factors:
      real*8, parameter ::
! eff is now a dbparam
!     &      eff=.7d0/3d0        ! locality+conversion efficiency of dissip -> kv
     &      n_column_ref=.01d0  ! col. N (1/s) at which kv z-scale = zscale_ref
     &     ,zscale_ref=100d0    ! reference vertical decay scale of kv (m)
     &     ,zscale_max=1000d0   ! maximum allowed vertical decay scale (m)
     &     ,delz_n_bot=300d0    ! distance over which near-bottom N is evaluated
     &     ,tpow_n_min=1d-5     ! lower bound on N (1/s) from offline dissipation
     &     ,rescale_tpow_max=2d0 ! maximum allowed rescaling factor for dissip rate
     &     ,dtnom=1d5           ! timestep (s) over which to apply trial kv
     &     ,kvnom=1d-5          ! peak magnitude of trial kv (m2/s)
     &     ,kv_max=10d0         ! maximum allowed output kv (m2/s)
     &     ,rescale_kv_max=kv_max/kvnom
c
      real*8 :: volgsp ! external func
      real*8 :: n_column,zscale,tpow_rescaled,n_bot,pbot,ztop,zbot,delz
      integer :: l,iter,ltop
      real*8 :: pesum0,pesum,rhokvbydze,vol,gl,sl,rescale_factor
      real*8, dimension(:), allocatable :: ze,fxs,fxg,rho

      allocate(ze(0:lm),fxg(0:lm),fxs(0:lm),rho(lm))

      kv(:) = 0.
      ! Two iterations.
      ! The first evaluates initial PE and calculates trial kv and its fluxes.
      ! The second evaluates final PE and rescales kv as per subr. header
      do iter=1,2
        ! compute potential energy
        pesum = 0.
        ze(lm) = 0.
        do l=lm,1,-1
          gl = g(l)
          sl = s(l)
          if(iter.eq.2) then ! apply flux convergence from trial kv
            gl = gl + dtnom*(fxg(l-1)-fxg(l))/m(l)
            sl = sl + dtnom*(fxs(l-1)-fxs(l))/m(l)
          endif
          vol = volgsp(gl,sl,p(l))
          ! subgrid quadrature is overkill to evaluate change of PE
          ze(l-1) = ze(l) + m(l)*vol
          pesum = pesum + .5*(ze(l-1)+ze(l))*m(l)
          rho(l) = 1d0/vol
        enddo
        pesum = pesum*grav
        if(iter.eq.1) then ! compute trial kv and evaluate its fluxes
          pesum0 = pesum ! save original column PE
          delz = .5*(ze(0)+ze(1)-ze(lm-1)-ze(lm))
          n_column = sqrt(grav*
     &         abs(volgsp(g(lm),s(lm),0d0)/
     &             volgsp(g(1),s(1),0d0)-1d0)/delz) + 1d-30
          zscale = min(zscale_ref*n_column_ref/n_column, zscale_max)
          do l=1,lm-1
            kv(l) = exp(-ze(l)/zscale)*kvnom
            rhokvbydze = kv(l)*(rho(l+1)+rho(l))/(ze(l-1)-ze(l+1))
            ! simplest possible diff. numerics should suffice here
            fxg(l) = (g(l)-g(l+1))*rhokvbydze
            fxs(l) = (s(l)-s(l+1))*rhokvbydze
          enddo
          fxg(0) = 0.
          fxs(0) = 0.
          fxg(lm) = 0.
          fxs(lm) = 0.
        elseif(pesum.eq.pesum0) then
c          write(6,*) 'pesum,pesum0 ',pesum,pesum0
c          write(6,*) 'ze ',ze
c          write(6,*) 'g ',g
c          write(6,*) 's ',s
c          write(6,*) 'm ',m
cc            call stop_model('testing',255)
          kv = 0.
        else
          ! rescale reference dissipation rate to get
          ! actual dissipation rate from near-bottom N
          zbot = .5d0*(ze(lm-1)+ze(lm))
          do ltop=lm-1,1,-1
            ztop = .5d0*(ze(ltop-1)+ze(ltop))
            delz = ztop-zbot
            if(ltop.eq.1 .or. delz .ge. delz_n_bot) exit
          enddo
          pbot = grav*sum(m)
          n_bot = sqrt(grav*
     &         abs(volgsp(g(lm),s(lm),pbot)/
     &             volgsp(g(ltop),s(ltop),pbot)-1d0) / delz )
          rescale_factor = n_bot/max(tpow_n_min,tpow_n)
          rescale_factor = min(rescale_tpow_max,rescale_factor)
          tpow_rescaled = tpow*rescale_factor
          ! (final kv) -> (trial kv) *
          !   (column dissip)/(rate of PE increase with trial kv)
          rescale_factor = dtnom*eff*tpow_rescaled/(pesum-pesum0)
          kv(:) = kv(:)*min(rescale_kv_max,max(0d0,rescale_factor))
          if(do_debug) then
            write(6,*) 'pesum,pesum0 debug',pesum,pesum0,
     &           (pesum-pesum0)/dtnom
            write(6,*) 'bv ',n_column
            write(6,*) 'kv ',kv
            write(6,*) 'tpow_n ',tpow_n,n_bot
          endif
        endif
      enddo
      deallocate(ze,fxg,fxs,rho)

      end subroutine get_kvtdiss

      SUBROUTINE STCONV
!@sum  STCONV uses vertical diffusion coefficients from KPP schmeme
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : grav,omega
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
#endif
      USE OCEAN,only : lmo,dts,ze,sinpo
      USE STRAITS, only : must,mmst,g0mst,gzmst,gxmst,s0mst,szmst,sxmst
     *     ,lmst,nmst,dist,wist,jst
#ifdef OCN_GISS_TURB
     *     ,otkest
#endif
#ifdef TRACERS_OCEAN
     *     ,trmst,txmst,tzmst
#endif
      USE ODIAG, only : olnst,ln_kvm,ln_kvg,ln_wgfl,ln_wsfl
#ifdef OCN_GISS_TURB
     &     ,ln_ri,ln_rrho,ln_bv2,ln_otke,ln_buoy
      USE GISSMIX_COM, only : otke
#endif
      IMPLICIT NONE
      REAL*8 MMLT,MML0
      REAL*8, DIMENSION(LMO,2) :: UL,G0ML,S0ML,GZML,SZML
      REAL*8, DIMENSION(LMO) :: MML,BYMML,DTBYDZ,BYDZ2,UL0,G0ML0,S0ML0
      REAL*8, DIMENSION(0:LMO+1) :: AKVM,AKVG,AKVS,AKVC
#ifdef OCN_GISS_TURB
     &     ,ri1,rrho,bv2,buoy
#endif
      REAL*8, DIMENSION(LMO) :: G,S,TO,BYRHO,RHO,PO,GHAT,FLG,FLS
      REAL*8, DIMENSION(LMO) :: DTP4
#ifdef TRACERS_OCEAN
      REAL*8 TRML(LMO,tracerlist%getsize(),2),
     &  TZML(LMO,tracerlist%getsize(),2),FLT(LMO,tracerlist%getsize())
      type(ocn_tracer_entry), pointer :: entry
      INTEGER ITR,NSIGT
#endif
C**** CONV parameters: BETA controls degree of convection (default 0.5).
      REAL*8, PARAMETER :: BETA=5d-1,BYBETA=1d0/BETA
#ifdef OCN_GISS_TURB
      LOGICAL, PARAMETER :: LDD = .true.
#else
      LOGICAL, PARAMETER :: LDD = .false.
#endif
      REAL*8, SAVE :: zgrid(0:LMO+1),hwide(0:LMO+1),byhwide(0:LMO+1)
      REAL*8 Shsq(LMO),dVsq(LMO)
     *     ,talpha(LMO),sbeta(LMO),galpha(LMO)
     &     ,dbloc(LMO),dbsfc(LMO),Ritop(LMO),
     *     alphaDT(LMO),betaDS(LMO),alphaDG(LMO)
      REAL*8, PARAMETER :: epsln=1d-20
      REAL*8, SAVE :: Bo,Bosol,bydts,ustar
      REAL*8 U2rho,RI,Coriol,HBL,HBLP,RHOM,RHO1,R,R2,DTBYDZ2,mld,rhoe
      REAL*8 VOLGSP,ALPHAGSP,BETAGSP,TEMGSP,SHCGS
      INTEGER, SAVE :: IFIRST = 1
      INTEGER I,L,N,LMIJ,IQ,ITER,NSIGG,NSIGS,KBL
#ifdef OCN_GISS_TURB
      REAL*8 rib(lmo),bf,ustarb2,exy
      INTEGER strait
c     REAL*8 wtnl(lmo)   !@var wtnl non-local term of wt
c     REAL*8 wsnl(lmo)   !@var wsnl non-local term of ws
      REAL*8 e(lmo)      !@var e ocean turbulent kinetic energy (m/s)**2
c     real*8 buoynl(lmo) !@var buoynl nonlocal part of buoy production
c     real*8 ga          !@var ga grav*alphaT
c     real*8 gb          !@var gb grav*alphaS
#endif

      IF (IFIRST.eq.1) THEN
        IFIRST=0
        zgrid(0) = epsln
        hwide(0) = epsln
        byhwide(0) = 0.
        DO L=1,LMO
          zgrid(L) = -0.5*(ZE(L-1) + ZE(L)) ! tracer level
          hwide(L) = zgrid(L-1) - zgrid(L) ! tracer level
          byhwide(L) = 1d0/hwide(L)
        END DO
        zgrid(LMO+1) = -ZE(LMO)
        hwide(LMO+1) = epsln
        byhwide(LMO+1) = 0.

        Ustar = 0.   ! No wind/ice stress in straits
        Bo = 0       ! No buoyancy forcing
        Bosol = 0
        BYDTS = 1d0/DTS
      END IF
C****
C**** Outside loop over straits
C****
      DO 790 N=1,NMST
C****
C**** Define parameters and currents
C****
      LMIJ=LMST(N)
      PO(1) = 5d-1*GRAV*MMST(1,N)/(DIST(N)*WIST(N))
      DO 220 L=1,LMIJ-1
      PO(L+1) = PO(L  ) + 5d-1*GRAV*(MMST(L,N)+MMST(L+1,N))
     *       /(DIST(N)*WIST(N))
      DTBYDZ(L) = DTS*DIST(N)*WIST(N)/MMST(L,N)
      BYDZ2(L)  = 2d0*DIST(N)*WIST(N)/(MMST(L+1,N)+MMST(L,N))
      MML(L)    = 5d-1*MMST(L,N)
  220 BYMML(L)  =  1d0/MML(L)
      MML(LMIJ) = 5d-1*MMST(LMIJ,N)
      BYMML(LMIJ) = 1d0/MML(LMIJ)
      DTBYDZ(LMIJ) = DTS*DIST(N)*WIST(N)/MMST(LMIJ,N)
C**** Loop over half boxes
      IQ=1
  240 RI = BETA*(2d0*IQ-3d0)
      DO L=1,LMIJ
        UL(L,IQ) = 5d-1*  MUST(L,N)*DIST(N)*BYMML(L)
      G0ML(L,IQ) = 5d-1*(G0MST(L,N) + RI*GXMST(L,N))
      S0ML(L,IQ) = 5d-1*(S0MST(L,N) + RI*SXMST(L,N))
      GZML(L,IQ) = 5d-1* GZMST(L,N)
      SZML(L,IQ) = 5d-1* SZMST(L,N)
#ifdef TRACERS_OCEAN
      TRML(L,:,IQ)=5d-1*(TRMST(L,N,:) + RI*TXMST(L,N,:))
      TZML(L,:,IQ)=5d-1* TZMST(L,N,:)
#endif
      END DO
C****
C**** Vertical mixing derived from KPP scheme
C****
  500 IF(LMIJ.LE.1)  GO TO 700

C**** Coriolis parameter, defined at strait end
      Coriol = 2d0*OMEGA*SINPO(JST(N,IQ))
C**** Save original values
        UL0 = UL(1:LMO,IQ)
      G0ML0 = G0ML(1:LMO,IQ)
      S0ML0 = S0ML(1:LMO,IQ)
C**** Start iteration for diffusivities
      ITER = 0
      HBL = 0
  510 ITER =ITER + 1
      HBLP = HBL
C**** Initiallise fields
      Shsq=0.
      dVsq=0.
      dbloc=0.
      dbsfc=0.
      Ritop=0.
      alphaDT=0. ; betaDS=0.

C**** velocity shears

      DO L=1,LMIJ-1
         Shsq(L)  = (UL(L,IQ)-UL(L+1,IQ))**2 !interface
         dvsq(L)  = (UL(1,IQ)-UL(L  ,IQ))**2 !tracer pts
      END DO
      dVsq(LMIJ) = (UL(1,IQ)-UL(LMIJ,IQ))**2 ! tracer pts

C****    density related quantities
C****
C****    density of surface layer                                (kg/m3)
C****          rho   = rho{G1),S(1),zt(1)}
C****    local buoyancy gradient at km interfaces:                (m/s2)
C****          dbloc = g/rho{k+1,k+1} * [ rho{k,k+1}-rho{k+1,k+1} ]
C****    buoyancy difference with respect to "zref", the surface  (m/s2)
C****          dbsfc = g/rho{k,k} * [ rho{k,k}-rho{1,k} ] CORRECT
C****    thermal expansion coefficient without 1/rho factor    (kg/m3/C)
C****          talpha= d(rho{k,k})/d(T(k))
C****    salt expansion coefficient without 1/rho factor     (kg/m3/PSU)
C****          sbeta =  d(rho{k,k})/d(S(k))

      G(1)     = G0ML(1,IQ)*BYMML(1)
      S(1)     = S0ML(1,IQ)*BYMML(1)
      BYRHO(1) = VOLGSP(G(1),S(1),PO(1))
      RHO(1)   = 1d0/BYRHO(1)
      DO L=2,LMIJ
         G(L)     = G0ML(L,IQ)*BYMML(L)
         S(L)     = S0ML(L,IQ)*BYMML(L)
         BYRHO(L) = VOLGSP(G(L),S(L),PO(L))
         RHO(L)   = 1d0/BYRHO(L)
         RHOM = 1d0/VOLGSP(G(L-1),S(L-1),PO(L)) ! L-1 density wrt PO(L)
         RHO1 = 1d0/VOLGSP(G(1),S(1),PO(L)) ! surf density wrt PO(L)
         dbsfc(L)   = GRAV*(1d0-RHO1*BYRHO(L))
         dbloc(L-1) = GRAV*(1d0-RHOM*BYRHO(L))
C**** numerator of bulk richardson number on grid levels
         Ritop(L) = (zgrid(1)-zgrid(L)) * dbsfc(L)
      END DO

C**** Double diffusive option (ddmix) needs alphaDT,betaDS
C**** alphaDT = mean -talpha * delta(temp.)    at interfaces  (kg/m3)
C**** betaDS  = mean sbeta  * delta(salt)     at interfaces  (kg/m3)

      if (LDD) then
        TO(1)      =    TEMGSP(G(1),S(1),PO(1))
        talpha(1) =  ALPHAGSP(G(1),S(1),PO(1)) ! <0
        sbeta(1)  =   BETAGSP(G(1),S(1),PO(1))*1d3
        do L=2,LMIJ
          TO(L)      =    TEMGSP(G(L),S(L),PO(L))
          talpha(L) =  ALPHAGSP(G(L),S(L),PO(L)) ! <0
          sbeta(L)  =   BETAGSP(G(L),S(L),PO(L))*1d3
          alphaDT(L-1) = - 5d-1 * (talpha(L-1) + talpha(L))
     $         * (TO(L-1) - TO(L))
          betaDS (L-1) = 5d-1 * (sbeta (L-1) + sbeta(L))
     $         * (S(L-1) - S(L))
        end do

        do l=1,lmij-1
          alphaDG(l) = rho(l+1) - 1d0/VOLGSP(G(L),S(L+1),PO(L+1))
          if(abs(g(l+1)-g(l)).gt.1d-10) then
            galpha(l) = alphaDG(l)/(g(l+1)-g(l))
          else
            galpha(l) = -1.
          endif
          betaDS(l) = rho(l+1) - 1d0/VOLGSP(G(L+1),S(L),PO(L+1))
          if(abs(s(l+1)-s(l)).gt.1d-10) then
            sbeta(l) = betaDS(l)/(s(l+1)-s(l))
          else
            sbeta(l) = 1000.
          endif
          rhoe = .5d0*(rho(l)+rho(l+1))
          galpha(l) = galpha(l)/rhoe
          sbeta(l) = sbeta(l)/rhoe
        enddo

      end if

#ifndef OCN_GISS_TURB
C**** Get diffusivities for the whole column
      CALL KPPMIX(LDD,ZE,zgrid,hwide,LMIJ,Shsq,dVsq,Ustar,Bo
     *     ,Bosol ,alphaDT,betaDS,dbloc,Ritop,Coriol,byhwide,
     *     AKVM,AKVS,AKVG,GHAT,HBL,KBL)
      akvc=akvs
#else
      call bldepth(ZE,zgrid,byhwide,LMIJ,dVsq,Ustar,Bo,Bosol
     *     ,dbloc,Ritop
     *     ,rib,HBL,KBL,bf)

      strait=1
      bf=0.d0
      ustarb2=0.d0
      exy=0.d0
      do l=1,lmij-1
         e(l)=otkest(l,n)
      end do

      call gissmix(
      ! in:
     &    lmij,ze,zgrid,dbloc,Shsq,alphaDT,betaDS,alphaDG,rho
     &   ,ustarb2,exy,Coriol,hbl,talpha,sbeta,galpha,strait,kbl,bf,ustar
     &   ,0,0
      ! out:
     &   ,ri1,rrho,bv2,akvm,akvg,akvs,akvc,e)
      buoy=0.
      do l=1,lmij-1
         otkest(l,n)=e(l)
         buoy(l)=-(akvg(l)-rrho(l)*akvs(l))/(1-rrho(l))*bv2(l)
      end do
#endif

C**** Correct units for diffusivities (m^2/s) => (kg^2/m^4 s)
      DO L=1,LMIJ-1
         R = 5d-1*(RHO(L)+RHO(L+1))
         R2 = R**2
         AKVM(L) = AKVM(L)*R2
         AKVG(L) = AKVG(L)*R2
         AKVS(L) = AKVS(L)*R2
         AKVC(L) = AKVC(L)*R2
         GHAT(L) = 0.  ! no non-local transports since no surface fluxes
         DTP4(L) = 0.
      END DO
      DTP4(LMIJ) = 0.

C**** For each field (U,G,S + TRACERS) call OVDIFF
C**** Momentum
      CALL OVDIFF(UL(1,IQ),AKVM(1),GHAT,DTP4
     &     ,ZE(1),ZGRID(1),DTBYDZ,BYDZ2,LMIJ,UL0)

C**** Enthalpy
      CALL OVDIFFS(G0ML(1,IQ),AKVG(1),GHAT,DTP4,DTBYDZ,BYDZ2,DTS
     *     ,LMIJ,G0ML0,FLG)
C**** Salinity
      CALL OVDIFFS(S0ML(1,IQ),AKVS(1),GHAT,DTP4,DTBYDZ,BYDZ2,DTS
     *     ,LMIJ,S0ML0,FLS)
      IF ((ITER.eq.1  .or. ABS(HBLP-HBL).gt.(ZE(KBL)-ZE(KBL-1))*0.25)
     *     .and. ITER.lt.4) GO TO 510
#ifdef TRACERS_OCEAN
C**** Tracers are diffused after iteration (GHAT always zero)
      DO ITR = 1,tracerlist%getsize()
        CALL OVDIFFS(TRML(1,ITR,IQ),AKVC(1),GHAT,DTP4,DTBYDZ,BYDZ2
     *       ,DTS,LMIJ,TRML(1,ITR,IQ),FLT(1,ITR))
      END DO
#endif
C**** Implicitly apply interpolated KV to linear profile
C**** No surface fluxes
      DTBYDZ2 = 12d0*DTBYDZ(1)**2*BYDTS
      GZML(1,IQ)=(GZML(1,IQ)+3d0*FLG(1))/(1d0+DTBYDZ2*AKVG(1))
      SZML(1,IQ)=(SZML(1,IQ)+3d0*FLS(1))/(1d0+DTBYDZ2*AKVS(1))
#ifdef TRACERS_OCEAN
      TZML(1,:,IQ)=(TZML(1,:,IQ)+3d0*FLT(1,:))/(1d0+DTBYDZ2*AKVC(1))
#endif
      DO L=2,LMIJ-1
        DTBYDZ2 = 6d0*DTBYDZ(L)**2*BYDTS
        GZML(L,IQ)=(GZML(L,IQ)+3d0*(FLG(L-1)+FLG(L)))
     *       /(1d0+DTBYDZ2*(AKVG(L-1)+AKVG(L)))
        SZML(L,IQ)=(SZML(L,IQ)+3d0*(FLS(L-1)+FLS(L)))
     *       /(1d0+DTBYDZ2*(AKVS(L-1)+AKVS(L)))
#ifdef TRACERS_OCEAN
        TZML(L,:,IQ)=(TZML(L,:,IQ)+3d0*(FLT(L-1,:)+FLT(L,:)))
     *       /(1d0+DTBYDZ2*(AKVC(L-1)+AKVC(L)))
#endif
      END DO
      DTBYDZ2 = 12d0*DTBYDZ(LMIJ)**2*BYDTS
      GZML(LMIJ,IQ)=(GZML(LMIJ,IQ)+3d0*FLG(LMIJ-1))
     *     /(1d0+DTBYDZ2*AKVG(LMIJ-1))
      SZML(LMIJ,IQ)=(SZML(LMIJ,IQ)+3d0*FLS(LMIJ-1))
     *     /(1d0+DTBYDZ2*AKVS(LMIJ-1))
#ifdef TRACERS_OCEAN
      TZML(LMIJ,:,IQ)=(TZML(LMIJ,:,IQ)+3d0*FLT(LMIJ-1,:))
     *     /(1d0+DTBYDZ2*AKVC(LMIJ-1))
#endif
C****
      DO L=1,LMIJ-1
        OLNST(L,N,LN_KVG) = OLNST(L,N,LN_KVG) + AKVG(L)
        OLNST(L,N,LN_KVM) = OLNST(L,N,LN_KVM) + AKVM(L)
#ifdef OCN_GISS_TURB
        OLNST(L,N,ln_ri)  = OLNST(L,N,ln_ri)  + ri1(L)
        OLNST(L,N,ln_rrho)= OLNST(L,N,ln_rrho)+ rrho(L)
        OLNST(L,N,ln_bv2) = OLNST(L,N,ln_bv2) + bv2(L)
        OLNST(L,N,ln_buoy)= OLNST(L,N,ln_buoy)+ buoy(L)
        OLNST(L,N,ln_otke)= OLNST(L,N,ln_otke)+ e(L)
#endif
C       OLNST(L,N,LN_KVS) = OLNST(L,N,LN_KVS) + AKVS(L)
        OLNST(L,N,LN_WGFL)= OLNST(L,N,LN_WGFL) + FLG(L)
        OLNST(L,N,LN_WSFL)= OLNST(L,N,LN_WSFL) + FLS(L)
      END DO

C**** End of half box loops
  700 IQ=IQ+1
      IF(IQ.LE.2)  GO TO 240
C****
C**** Recreate main prognostic variables of potential enthalpy and
C**** salt from those on the half boxes
C****
      DO L=1,LMIJ
       MUST(L,N) = (  UL(L,2) +   UL(L,1))*MML(L)/DIST(N)
      G0MST(L,N) =  G0ML(L,2) + G0ML(L,1)
      NSIGG = EXPONENT(G0MST(L,N)) -1 - 42
      GXMST(L,N) = (G0ML(L,2) - G0ML(L,1))*BYBETA
      CALL REDUCE_FIG(NSIGG,GXMST(L,N))
      GZMST(L,N) =  GZML(L,2) + GZML(L,1)
      S0MST(L,N) =  S0ML(L,2) + S0ML(L,1)
      NSIGS = EXPONENT(S0MST(L,N)) -1 - 42 + 4
      SXMST(L,N) = (S0ML(L,2) - S0ML(L,1))*BYBETA
      CALL REDUCE_FIG(NSIGS,SXMST(L,N))
      SZMST(L,N) =  SZML(L,2) + SZML(L,1)
C****  limit salinity gradients
      if ( abs(SXMST(L,N)) > S0MST(L,N) )
     *     SXMST(L,N) = sign(S0MST(L,N),SXMST(L,N))
      if ( abs(SZMST(L,N)) > S0MST(L,N) )
     *     SZMST(L,N) = sign(S0MST(L,N),SZMST(L,N))
C****
#ifdef TRACERS_OCEAN
      DO ITR=1,tracerlist%getsize()
        entry=>tracerlist%at(itr)
        TRMST(L,N,ITR) = TRML(L,ITR,2) + TRML(L,ITR,1)
        NSIGT = EXPONENT(TRMST(L,N,ITR)) -1 - 42
        TXMST(L,N,ITR) =(TRML(L,ITR,2) - TRML(L,ITR,1))*BYBETA
        CALL REDUCE_FIG(NSIGT,TXMST(L,N,ITR))
        TZMST(L,N,ITR) = TZML(L,ITR,2) + TZML(L,ITR,1)
C****
        if (entry%t_qlimit) then  ! limit gradients
          if ( abs(TXMST(L,N,ITR)) > TRMST(L,N,ITR) )
     *         TXMST(L,N,ITR) = sign(TRMST(L,N,ITR),TXMST(L,N,ITR))
          if ( abs(TZMST(L,N,ITR)) > TRMST(L,N,ITR) )
     *         TZMST(L,N,ITR) = sign(TRMST(L,N,ITR),TZMST(L,N,ITR))
        end if
C****
      END DO
#endif
      END DO
C**** End of outside loop over straits
  790 CONTINUE
      RETURN
      END SUBROUTINE STCONV


      SUBROUTINE OVDIFF(U,K,GHAT,DTP4,ZE,Z,DTBYDZ,BYDZ2,LMIJ,U0)
!@sum  OVDIFF Implicit vertical diff + non local transport for velocity
!@auth Gavin Schmidt
      USE OCEAN, only : LMO
      USE TRIDIAG_MOD, only : tridiag
      IMPLICIT NONE
      REAL*8, DIMENSION(LMO), INTENT(IN) :: U0,K,GHAT,DTBYDZ,BYDZ2
     &                                     ,ZE,Z,DTP4 ! Z<0
      REAL*8, DIMENSION(LMO), INTENT(OUT) :: U
      INTEGER, INTENT(IN) :: LMIJ
      REAL*8, DIMENSION(LMO) :: A,B,C,R
      INTEGER L
      REAL*8 tmp
C****  U0,U  input and output field (velocity or concentration)
C****     K  vertical diffusivity ((z)^2/ s)
C****  GHAT  non local transport of scalar
C****           kv * ghats * surface flux
C**** DTBYDZ  DT/DZ_L
C**** BYDZ2  1d0/DZ_L+1/2
C****    DT  timestep (s)
C**** top boundary::
C**** U(1)-U0(1)=-DTBYDZ(1)*(BYDZ2(1)*K(1)*(U(1)-U(2))-GHAT(1))+DTP4   
C**** Boundary conditions assumed to be no-flux at Z=0, Z=Z(LMIJ)
C**** Calculate operators for tridiagonal solver
      A(1) = 0
c     B(1) = 1d0   + DTBYDZ(1)*BYDZ2(1)*K(1)
c     C(1) =       - DTBYDZ(1)*BYDZ2(1)*K(1)
      C(1) =       - DTBYDZ(1)*BYDZ2(1)*K(1)
      B(1) = 1d0-C(1)   
      
      tmp  = U0(1)
c     tmp  = .5d0*(U0(1)+U0(2))
      R(1) = tmp - DTBYDZ(1)*GHAT(1) + DTP4(1)
c     R(1) = tmp
      DO L=2,LMIJ-1
c       A(L) =       - DTBYDZ(L)* BYDZ2(L-1)*K(L-1)
        A(L)=-DTBYDZ(L)*BYDZ2(L-1)*K(L-1)
        C(L)=-DTBYDZ(L)*BYDZ2(L)*K(L)
c       A(L)=-DTBYDZ(L)*BYDZ2(L-1)*K(L-1)
c       C(L)=-DTBYDZ(L)*BYDZ2(L)*K(L)
c       B(L)=1-(A(L)+C(L))
        B(L) = 1d0   + DTBYDZ(L)*(BYDZ2(L-1)*K(L-1)+BYDZ2(L)*K(L))
c       C(L) =       - DTBYDZ(L)*                   BYDZ2(L)*K(L)
        tmp  = U0(L)
c       tmp  = .5d0*(U0(L-1)+U0(L+1))
        R(L) = tmp + DTBYDZ(L)*(GHAT(L-1) - GHAT(L)) + DTP4(L)
      END DO
      A(LMIJ) =          - DTBYDZ(LMIJ)*BYDZ2(LMIJ-1)*K(LMIJ-1)
      B(LMIJ) = 1d0-A(LMIJ)      
      C(LMIJ) = 0
      tmp  = U0(LMIJ)
c     tmp  = .5d0*(U0(LMIJ-1)+U0(LMIJ))
      R(LMIJ) = tmp + DTBYDZ(LMIJ)*GHAT(LMIJ-1) + DTP4(LMIJ)
c     R(LMIJ) = tmp

      CALL TRIDIAG(A,B,C,R,U,LMIJ)

      RETURN
      END SUBROUTINE OVDIFF

      SUBROUTINE OVDIFFS(U,K,GHAT,DTP4,DTBYDZ,BYDZ2,DT,LMIJ,U0,FL)
!@sum  OVDIFFS Implicit vertical diff + non local transport for tracers
!@auth Gavin Schmidt
      USE OCEAN, only : LMO
      USE TRIDIAG_MOD, only : tridiag
      IMPLICIT NONE
      REAL*8, DIMENSION(LMO), INTENT(IN) :: U0,K,GHAT,DTP4,DTBYDZ,BYDZ2
      REAL*8, DIMENSION(LMO), INTENT(OUT) :: U,FL
      REAL*8, INTENT(IN) :: DT
      INTEGER, INTENT(IN) :: LMIJ
      REAL*8, DIMENSION(LMO) :: A,B,C,R
      INTEGER L
C****  U0,U  input and output fields (total tracer)
C****     K  vertical diffusivity ((z)^2/ s)
C****  GHAT  non local transport of scalar
C****            kv * ghats * surface flux
C**** DTBYDZ  DT/DZ_L
C**** BYDZ2  1d0/DZ_L+1/2
C****    DT  timestep (s)
C****    FL  diffusive flux at boundary in units of total tracer
C**** Boundary conditions assumed to be no-flux at Z=0, Z=Z(LMIJ)
C**** Calculate operators for tridiagonal solver
      A(1) = 0
      B(1) = 1d0   + DTBYDZ(1)*BYDZ2(1)*K(1)
      C(1) =       - DTBYDZ(2)*BYDZ2(1)*K(1)
      R(1) = U0(1) - DT * GHAT(1) + DTP4(1)
      DO L=2,LMIJ-1
        A(L) =       - DTBYDZ(L-1)* BYDZ2(L-1)*K(L-1)
        B(L) = 1d0   + DTBYDZ(L  )*(BYDZ2(L-1)*K(L-1)+BYDZ2(L)*K(L))
        C(L) =       - DTBYDZ(L+1)*                   BYDZ2(L)*K(L)
        R(L) = U0(L) + DT * (GHAT(L-1) - GHAT(L)) + DTP4(L)
      END DO
      A(LMIJ) =          - DTBYDZ(LMIJ-1)*BYDZ2(LMIJ-1)*K(LMIJ-1)
      B(LMIJ) = 1d0      + DTBYDZ(LMIJ  )*BYDZ2(LMIJ-1)*K(LMIJ-1)
      C(LMIJ) = 0
      R(LMIJ) = U0(LMIJ) + DT * GHAT(LMIJ-1) + DTP4(LMIJ-1)

      CALL TRIDIAG(A,B,C,R,U,LMIJ)

      DO L=1,LMIJ-1
        FL(L)=K(L)*(DTBYDZ(L+1)*U(L+1)-DTBYDZ(L)*U(L))*BYDZ2(L)
     &       -DT*GHAT(L) ! include nonlocal part
      END DO
C****
      RETURN
      END SUBROUTINE OVDIFFS

      SUBROUTINE REDUCE_FIG(NSIG,RX)
!@sub reduce_fig reduce significant figures if calculation is garbage
!@+   made separate to avoid OMP compiler bug (due to NINT)
!@auth Gavin Schmidt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NSIG
      REAL*8, INTENT(INOUT) :: RX

      IF (NSIG+30.gt.EXPONENT(RX)) RX =
     *     SCALE(REAL(NINT(SCALE(RX,-NSIG)),KIND=8),NSIG)

      RETURN
      END SUBROUTINE REDUCE_FIG

      SUBROUTINE alloc_kpp_com(grid)
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Reto Ruedy
      use pario, only : par_open,read_dist_data,par_close
      USE DOMAIN_DECOMP_1D, only : dist_grid,getDomainBounds
!      USE OCEANR_DIM

      USE KPP_COM
      use dictionary_mod, only : sync_param
#ifdef TRACERS_OCEAN
      use ocn_tracer_com, only: tracerlist
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER
      integer :: fid

      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

c     ALLOCATE(  KPL(IM,J_0H:J_1H)    , STAT = IER)
#ifdef OCN_GISS_SM
c     ALLOCATE(  ktap(IM,J_0H:J_1H)    , STAT = IER)
c     ALLOCATE(  k02count(IM,J_0H:J_1H), STAT = IER)
c     ktap=1
c     k02count=1.
#endif

      ALLOCATE( G0M1(IM,J_0H:J_1H,LSRPD) , STAT = IER)

      ALLOCATE(  MO1(IM,J_0H:J_1H),    S0M1(IM,J_0H:J_1H),
     *          GXM1(IM,J_0H:J_1H),    SXM1(IM,J_0H:J_1H),
     *          GYM1(IM,J_0H:J_1H),    SYM1(IM,J_0H:J_1H),
     *           UO1(IM,J_0H:J_1H),     VO1(IM,J_0H:J_1H),
     *           UOD1(IM,J_0H:J_1H),    VOD1(IM,J_0H:J_1H),
     *   STAT = IER)
#ifdef TRACERS_OCEAN
      ALLOCATE( TRMO1(tracerlist%getsize(),IM,J_0H:J_1H),
     *          TXMO1(tracerlist%getsize(),IM,J_0H:J_1H),
     *          TYMO1(tracerlist%getsize(),IM,J_0H:J_1H),
     *   STAT = IER)
#endif

      call sync_param("ocean_use_tdiss",use_tdiss)
      if(use_tdiss==1) then
#ifdef ENHANCED_DEEP_MIXING
        call stop_model(
     &    'ENHANCED_DEEP_MIXING and use_tdiss==1 incompatible',255)
#endif
        allocate(tdiss(im,j_0h:j_1h), stat = ier)
        fid = par_open(grid,'TDISS','read')
        call read_dist_data(grid,fid,'dissip',tdiss)
        call par_close(grid,fid)
        allocate(tdiss_n(im,j_0h:j_1h), stat = ier)
        fid = par_open(grid,'TDISS_N','read')
        call read_dist_data(grid,fid,'buoyancy',tdiss_n)
        call par_close(grid,fid)

        call sync_param("ocean_tdiss_eff",tdiss_eff)
      endif

      END SUBROUTINE alloc_kpp_com

      subroutine get_gradients0(mokgm2,q_in,flag,qx,qy)
      use domain_decomp_1d, only :
     &     getDomainBounds,halo_update,south,north
      use oceanr_dim, only : grid=>ogrid
      use ocean, only : dxpo,dyvo,dxypo
      use ocean, only : lmu,lmm,
     &     nbyzm,nbyzu,nbyzv, i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      use ocean, only : im,jm,lmo,ivnp,sinic,cosic,sinu,cosu
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     mokgm2, ! units: kg/m2
     &     q_in, ! extensive or intesive units
     &     qx,qy ! outputs have intensive units
      integer flag ! 1: q_in is extensive; 0: q_in is intensive

      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     q,qxe,qyn
      integer :: i,j,l,n
      integer :: j_0s,j_1s,j_0,j_1
      logical :: have_north_pole
      real*8 :: unp,vnp

      call getdomainbounds(grid,
     &     j_strt=j_0, j_stop=j_1,
     &     j_strt_skp=j_0s, j_stop_skp=j_1s,
     &     have_north_pole=have_north_pole)

      if(flag.eq.1) then ! q_in is extensive
        ! convert q to intensive units
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          q(i,j,l) = q_in(i,j,l)/(mokgm2(i,j,l)*dxypo(j))
        enddo
        enddo
        enddo
        enddo
      else               ! q_in is intensive
        q=q_in
      endif

      call halo_update(grid,q,from=north)

      if(have_north_pole) then
        do l=1,lmo
          q(2:im,jm,l) = q(1,jm,l)
        enddo
      endif

      ! compute gradients at cell edges.  gradients at coastlines
      ! are zero
      qxe = 0.
      qyn = 0.
      do l=1,lmo
        do j=j_0s,j_1s
          do n=1,nbyzu(j,l)
          do i=i1yzu(n,j,l),min(im-1,i2yzu(n,j,l))
            qxe(i,j,l) = (q(i+1,j,l)-q(i,j,l))/dxpo(j)
          enddo
          enddo
          i=im
          if(lmu(i,j).ge.l) then
            qxe(i,j,l) = (q(1,j,l)-q(i,j,l))/dxpo(j)
          endif
          do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            qyn(i,j,l) = (q(i,j+1,l)-q(i,j,l))/dyvo(j)
          enddo
          enddo
        enddo
      enddo

      ! average cell-edge gradients to cell centers.
      ! gradients in the north polar cap temporarily left at zero.
      call halo_update(grid,qyn,from=south)

      qx = 0.
      qy = 0.
      do l=1,lmo
      do j=j_0s,j_1s
        do n=1,nbyzm(j,l)
        do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
          qx(i,j,l) = .5*(qxe(i-1,j,l)+qxe(i,j,l))
          qy(i,j,l) = .5*(qyn(i,j-1,l)+qyn(i,j,l))
        enddo
        enddo
        i=1
        if(lmm(i,j).ge.l) then
          qx(i,j,l) = .5*(qxe(im,j,l)+qxe(i,j,l))
          qy(i,j,l) = .5*(qyn(i,j-1,l)+qyn(i,j,l))
        endif
      enddo
      enddo

      if(have_north_pole) then
c at the north pole
        unp = 0.
        vnp = 0.
        j = jm-1
        do l=1,lmo
          do n=1,nbyzv(j,l)
            do i=i1yzv(n,j,l),i2yzv(n,j,l)
              unp = unp - sinic(i)*qyn(i,j,l)
              vnp = vnp + cosic(i)*qyn(i,j,l)
            enddo
          enddo
          unp = unp*2/im
          vnp = vnp*2/im
          do i=1,im
            qx(i,jm,l) = unp*cosu(i)  + vnp*sinu(i)
            qy(i,jm,l) = vnp*cosic(i) - unp*sinic(i)
          enddo
c         qx(im,jm,l) = unp   ! as a result of the above loop
c         qx(ivnp,jm,l) = vnp ! as a result of the above loop
        enddo
      endif
      return
      end subroutine get_gradients0
