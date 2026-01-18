C****   
C**** OCEAN_COM.f    Module Variables for Ocean    2015/02/18
C****
#include "rundeck_opts.h"

      Module OCEAN
!@sum  OCEAN dynamic ocean related variables
!@auth Gary Russell/Gavin Schmidt
C**** Note that we currently use the same horizontal grid as for the
C**** atmosphere. However, we can redefine im,jm if necessary.
      Use CONSTANT,  Only: TWOPI
      Use OCEANRES,  Only: IM=>IMO,JM=>JMO, LMO, LMO_MIN, dZO
      Use SparseCommunicator_mod
#ifdef CUBED_SPHERE
      use cs2ll_utils, only : aoremap_type=>xgridremap_type
#else
      use hntrp_mod, only :   aoremap_type=>hntrp_type
#endif
      Implicit None
      Integer*4,Parameter ::
     *  IVSP = 3*IM/4,      !  V at south pole is stored in U(IVSP,1)
     *  IVNP =   IM/4       !  V at north pole is stored in U(IVNP,JM)
      REAL*8, PARAMETER :: FIM=IM, BYIM=1./FIM

!@dbparam OBottom_drag use ocean bottom drag routine (default=1)
!@dbparam OCoastal_drag use ocean coastal drag routine (default=1)
!@dbparam OTIDE: Lunar,Solar tides accelerate UO,VO,UOD,VOD (default=0)
      Integer*4 ::
     *     OBottom_drag = 1,    
     *     OCoastal_drag = 1,
     *     OTIDE = 0

!@dbparam oc_mean_salinity define mean salinity of ocean (if set)
!@dbparam oc_tracer_mean mean tracer ratio of ocean 
!@+       (in permil for water isotopes)
      REAL*8 :: 
     *     oc_salt_mean = -999. 
#ifdef TRACERS_OCEAN
      ! move to ocn_tracer_com?
      REAL*8, DIMENSION(:), ALLOCATABLE :: oc_tracer_mean
#endif 
!@var MO mass of ocean (kg/m^2)
!@var UO E-W velocity on C-grid (m/s)
!@var VO N-S velocity on C-grid (m/s)
!@var UOD E-W velocity on D-grid (m/s)
!@var VOD N-S velocity on D-grid (m/s)
!@var G0M,GXMO,GYMO,GZMO pot. enthalpy of ocean (+moments) (J)
!@var S0M,SXMO,SYMO,SZMO salinity of ocean (+moments) (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: MO,UO,VO,UOD,VOD,
     *     G0M,S0M

      INTEGER :: USE_QUS=0
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     GXMO,GYMO,GZMO, GXXMO,GYYMO,GZZMO, GXYMO,GYZMO,GZXMO,
     *     SXMO,SYMO,SZMO, SXXMO,SYYMO,SZZMO, SXYMO,SYZMO,SZXMO
#ifdef OCN_GISS_SM
     *    ,rx,ry,gx,gy,sx,sy
#endif

      Real*8 UONP(LMO), !  U component at north pole, points down 90W
     *       VONP(LMO)  !  V component at north pole, points down 0 (GM)

C**** ocean geometry (should this be in a separate module?)
      Real*8 
     *     DLON,    !@var DLON longitudinal spacing in radians
     *     DLAT,    !@var DLAT latitudinal spacing in radians
     *     DLATM,   !@var DLATM latitudinal spacing in minutes
     *     FJEQ     !@var FJEQ location of equator in grid units 
     *   , oDLAT_DG ! grid spacing in latitude (deg) 
     *   , oDLON_DG ! grid spacing in longitude (deg) 

      REAL*8, ALLOCATABLE, DIMENSION(:,:):: OXYP

      REAL*8, DIMENSION(JM) :: DXYPO,DXPO,DYPO,DXVO,DYVO
     *     ,COSPO,SINPO,DXYVO,DXYSO,DXYNO,RAMVS,RAMVN,RLAT,BYDXYPO
      REAL*8, DIMENSION(0:JM) :: SINVO,COSVO
      REAL*8, DIMENSION(0:LMO) :: ZE,DZOE,BYDZOE
      REAL*8, DIMENSION(LMO) :: ZMID
      Real*8
     *  DXPGF(0:JM),! DXYV/dYPGF is north-south distance used in PGF
     *  DYPGF(JM), !  DXYP/dXPGF is east-west distance used in PGF
     *  COSM(JM),  !  = .5*COSV(J-1) + .5*COSV(J)
     *  COSQ(JM),  !  sQuare of COSine = .5*COSV(J-1)^2 + .5*COSV(J)^2
     *  SINIC(IM), !  SINe of longitude of grid cell center from IDL
     *  COSIC(IM), !  COSine of longitude of grid cell center from IDL
     *  SINU(IM),  !  SINe of longitude of eastern edge of grid cell
     *  COSU(IM),  !  COSine of longitude of eastern edge of grid cell
     *  SINxY(JM), !  DLAT * SINe of latitude used by Coriolis force
     *  TANxY(JM)  !  DLAT * TANgent of latitude used by metric term
!@var  oLAT_DG latitude of mid points of primary and sec. grid boxs (deg)
      REAL*8, DIMENSION(JM,2) :: oLAT_DG
!@var  oLON_DG longitude of mid points of prim. and sec. grid boxes (deg)
      REAL*8, DIMENSION(IM,2) :: oLON_DG
!@var  IMAXJ varying number of used longitudes
      INTEGER, DIMENSION(JM) :: IMAXJ
      INTEGER, DIMENSION(:,:), allocatable :: LMM,LMU,LMV
      Integer*4 J40S, !  maximum grid cell below 40S (used in OPFIL)
     *           J1O  !  most southern latitude (J) where ocean exists
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oLAT2D_DG !distributed latitute array (in degrees)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: HATMO,HOCEAN,FOCEAN
C**** ocean related parameters
      INTEGER NDYNO,MDYNO,MSGSO
!@dbparam DTO timestep for ocean dynamics (s)
      REAL*8 :: DTO=450.        ! default. setable parameter
      REAL*8 DTOFS,DTOLF,DTS,BYDTS

!@dbparam NOCEAN number of ocean advective timesteps per physics timestep
!  NDYNO must be multiple of 2*NOCEAN
      integer :: NOCEAN = 1
     &     + JM/180   ! force default of 2 for 1-degree res

      ! to-do: make binomial filter coeffs into rundeck parameters
      ! if/when it is reinstated
      integer, parameter ::
     &      NORDER=4      !  order of Alternating Binomial Filter (must be even)
      ! coeffs for divergence/vorticity filter in X/Y dirs. 
      real*8, parameter ::
     &  OABFUX = 0d0
     &    + (jm/180)*(.15d0/4**norder),  ! switch on for 1-degree res
     &  OABFVX = OABFUX,
     &  by4tonv = oabfux,
     &  by4tonu = by4tonv

!@var budget grid quantities (defined locally on each proc.)
      REAL*8, ALLOCATABLE, DIMENSION(:,:):: owtbudg
      INTEGER, ALLOCATABLE, DIMENSION(:,:):: oJ_BUDG
      INTEGER :: oJ_0B,oJ_1B
!@var OPRESS Anomalous pressure at surface of ocean (under ice) (Pa)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OPRESS ! (IM,JM)
!@var OGEOZ ocean geopotential at surface (m^2/s^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OGEOZ,OGEOZ_SV
!@var KPL level to which mixed layer descends (1)
      integer, ALLOCATABLE, DIMENSION(:,:) :: kpl
!@var OPBOT ocean bottom pressure (diagnostic only) (Pa)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OPBOT

#ifdef TRACERS_OCEAN
!@var TRMO,TXMO,TYMO,TZMO tracer amount (+moments) in ocean (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRMO

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) ::
     *     TXMO,TYMO,TZMO, TXXMO,TYYMO,TZZMO, TXYMO,TYZMO,TZXMO
#endif
     
      type (SparseCommunicator_type), save :: mySparseComm_type

!@var nbyz[muvc]: # of basins at each lat/depth
!@var i[12]yz[muvc]: start/end i-indices for each basin
! m: cell center ! u: east edge ! v: north edge ! c: northeast corner
      integer, parameter :: nbyzmax=22 ! suffices up to 1x1.25 deg res
      integer, dimension(:,:), allocatable :: nbyzm,nbyzu,nbyzv,nbyzc
      integer, dimension(:,:,:), allocatable ::
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv, i1yzc,i2yzc

!@var remap_a2o,remap_o2a atm->ocn,ocn->atm interpolation info
      type(aoremap_type) ::
     &     remap_a2o            ! atm A -> ocn A
     &    ,remap_o2a            ! ocn A -> atm A

      REAL*8, SAVE, ALLOCATABLE, DIMENSION(:) ::
     *     BYDXYV, KHP,KHV,TANP,TANV,BYDXV,BYDXP,BYDYV,BYDYP
      REAL*8, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::
     *      UXA,UXB,UXC,UYA,UYB,UYC,VXA,VXB,VXC,VYA,VYB,VYC
      REAL*8, SAVE :: BYDXYPJM
      REAL*8, SAVE, DIMENSION(LMO) :: UYPB
      REAL*8, SAVE, DIMENSION(IM,LMO) :: UYPA
      REAL*8, PARAMETER :: FSLIP=0.

!@param zmax_glmelt max. nominal depth over which to deposit glacial melt (m)
! to-do: make this a rundeck parameter
      real*8, parameter :: zmax_glmelt = 202d0 ! taken from 32-layer version

! Various derived quantities for the convenience of routines needing them.
!@var g3d potential enthalpy (J/kg)
!@var t3d in-situ temperature (degC, ref to mid point pressure)
!@var s3d salinity (kg/kg)
!@var p3d mid point pressure (Pa)
!@var r3d in-situ density (kg/m3, ref to mid point pressure)
!@var v3d specific volume (m3/kg, ref to mid point pressure)
      real*8, allocatable, dimension(:,:,:) ::
     &     g3d,t3d,s3d,p3d,r3d,v3d

#ifdef TRACERS_OCEAN
!@dbparam ntrtrans number of physics timesteps per tracer transport
!@+   timestep.   Here, "tracer transport timestep" refers only to
!@+   advective and mesoscale transport.   Other transport proceses,
!@+   such as sedimentation and turbulent vertical mixing, are applied
!@+   to tracers (the trmo array et al.) every physics timestep regardless
!@+   of ntrtrans, as are "source terms" such as surface fluxes and ice
!@+   formation. See remarks below on the motr variable.
!@+   Current restrictions on this speedup method:
!@+     - mod(nday,ntrtrans) = mod(nssw,ntrtrans) = 0 for convenience.
!@+     - In E2.1, ntrtrans*dt > 3 hours occasionally crashes, but
!@+       there is currently no speedup benefit going past 3 hours.
!@+     - It should not be used for water tracers if those tracers
!@+       must remain exactly consistent with salinity, since salt
!@+       is still transported on physics timesteps.   Moreover,
!@+       some work remains to be done on isotope fractionation terms
!@+       in ground_oc.
      integer :: ntrtrans=1
!@var do_tracer_trans whether the current physics timestep is a
!@+   tracer transport timestep, ie whether mod(step,ntrtrans)=0
      logical :: do_tracer_trans
!
!@var asm[uvw] advective mass fluxes accumulated over the ntrtrans
!@+   physics timesteps in a tracer transport timestep.
!
!@var motr (kg/m2/layer) the current seawater mass corresponding to tracer
!@+   state.
!@+   motr is equal to mo immediately after tracer transports are
!@+   performed every ntrtrans physics timesteps, but motr will then
!@+   differ from mo until the next transport timestep, since source
!@+   terms are the only tendencies applied to motr during the interim,
!@+   while mo evolves via all tendencies each physics timestep
!@+   (including advective mass fluxes and bolus terms).
!@+   Schematic of the cycle:
!@+     (1) for each of ntrtrans physics timesteps
!@+         (a) mo,G,S = mo,G,S + all increments for this step
!@+         (b) motr,trmo = motr,trmo + all increments for this step
!@+                    excluding mesoscales and advection
!@+         (c) accumulate mass fluxes to be applied in (2)
!@+     (2) perform long-step transport using 1c. motr differs from
!@+         mo before this operation and is exactly equal to mo afterward.
!@+     (3) back to (1) for next cycle
!@+   To diagnose the mixing ratio q of a tracer during the steps (1),
!@+   motr rather than mo must be used as the denominator in
!@+   q = (tracer mass)/(seawater mass), ie q = trmo/motr.
!@+   With ntrtrans not more than a few hours, instantaneous differences
!@+   between mo and motr during (1) are assumed to be small enough that
!@+   quantities defined w.r.t. mo, such as temperature, salinity, and
!@+   turbulent vertical diffusivity, can be applied to the slightly
!@+   different volumes occupied by motr.
!@+   TODO: make some decisions about how to treat mo versus motr in the
!@+   monthly output diagnostics (e.g. mo->motr in parts of denom_toijl).
!
!@var mosv0 (kg/m2/layer) a book-keeping array to facilitate adding source
!@+   terms to motr
      real*8, allocatable, dimension(:,:,:) :: asmu,asmv,asmw,motr,mosv0
#endif

      contains

      subroutine alloc_odiff(grid)
      use DOMAIN_DECOMP_1D, only: dist_grid, getDomainBounds
      type (dist_grid) :: grid

      integer :: J_0H, J_1H

      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      allocate( BYDXYV(grid%j_strt_halo:grid%j_stop_halo) )
      allocate( KHP   (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( KHV   (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( TANP  (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( TANV  (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BYDXV (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BYDXP (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BYDYV (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BYDYP (grid%j_strt_halo:grid%j_stop_halo) )
      
      allocate( UXA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( UXB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( UXC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( UYA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( UYB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( UYC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VXA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VXB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VXC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VYA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VYB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VYC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

      end subroutine alloc_odiff

      END Module OCEAN

      Module OCEAN_DYN
!@sum  OCEAN_DYN contains variables used in ocean dynamics
!@auth Gavin Schmidt/Gary Russell
      Use OCEAN, Only : im,jm,lmo
!@var DH height of each ocean layer
!@var VBAR mean specific volume of each layer
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DH,VBAR !  (IM,JM,LMO)
     &     ,dZGdP,BYDH

!@var GUP,GDN specific pot enthropy upper,lower part of layer (J/kg)
!@var SUP,SDN salinity at           upper,lower part of layer (1)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GUP,GDN,SUP,SDN

C**** momentum and mass fluxes
!@var MMI initial mass field (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: MMI  !  (IM,JM,LMO)
!@var SMU,SMV,SMW integrated mass fluxes
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SMU,SMV,SMW ! (IM,JM,LMO)
!@var CONV mass flux convergence
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CONV      !   (IM,JM,LMO)
!@var MU,MV,MW instantaneous mass fluxes
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: MU,MV,MW  !   (IM,JM,LMO)
C****
      END Module OCEAN_DYN

      Module SW2OCEAN
!@sum  SW2OCEAN variables for putting solar radiation into ocean
!@auth Gavin Schmidt/Gary Russell
      Use OCEAN, Only : ze,lmo
      Implicit None

!@param zmax_solar nominal maximum SW penetration depth (m)
      real*8, parameter :: zmax_solar=92d0 ! taken from 32-layer version
      integer :: lsrpd=0 ! layer (bottom) index corresponding to zmax_solar

      real*8, dimension(:), allocatable :: FSR,FSRZ,dFSRdZ,dFSRdZB

      REAL*8, PARAMETER :: RFRAC=.62d0, ZETA1=1.5d0, ZETA2=2d1

      CONTAINS

      SUBROUTINE init_solar
!@sum  init_solar calculates penetration of solar radiation
!@auth Gavin Schmidt/Gary Russell
      REAL*8 EF,EFZ,Z
      INTEGER L
      EF(Z) = RFRAC*EXP(-Z/ZETA1) + (1d0-RFRAC)*EXP(-Z/ZETA2)
      EFZ(Z)=ZETA1*RFRAC*EXP(-Z/ZETA1)+ZETA2*(1d0-RFRAC)*EXP(-Z/ZETA2)

      ! determine lsrpd from zmax_solar and layering
      do lsrpd=1,lmo-1
        if(ze(lsrpd+1) > zmax_solar) exit
      enddo
      allocate(fsr(lsrpd),fsrz(lsrpd),dfsrdz(lsrpd),dfsrdzb(lsrpd))

C****
C**** Calculate the fraction of solar energy absorbed in each layer
C****
      do l=1,LSRPD
         FSR(l) = EF(ZE(l-1))
         FSRZ(l) = -3d0*(EF(ZE(l-1))+EF(ZE(l))) +
     *        6d0*(EFZ(ZE(l-1))-EFZ(ZE(l)))/(ZE(l)-ZE(l-1))
         dFSRdZB(l) = FSR(l)/(ZE(l)-ZE(l-1))
      end do
      do l=1,LSRPD-1
         dFSRdZ(l) = (FSR(l)-FSR(l+1))/(ZE(l)-ZE(l-1))
      end do
      dFSRdZ (LSRPD) = 0.
C****
      END SUBROUTINE init_solar
C****
      END Module SW2OCEAN

      SUBROUTINE alloc_ocean
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Rodger Abel
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, am_i_root
      USE OCEANR_DIM, only : ogrid,J_0H,J_1H,init_oceanr_grid  

      USE OCEANRES, only : IM=>IMO, JM=>JMO, LMO 

      USE OCEAN, only : MO,G0M,S0M
      USE OCEAN, only : UO,VO,UOD,VOD
      USE OCEAN, only : OPRESS,OPBOT, OGEOZ,OGEOZ_SV,kpl
      USE OCEAN, only : use_qus,
     *     GXMO,GYMO,GZMO, GXXMO,GYYMO,GZZMO, GXYMO,GYZMO,GZXMO,
     *     SXMO,SYMO,SZMO, SXXMO,SYYMO,SZZMO, SXYMO,SYZMO,SZXMO,
     &     g3d,t3d,s3d,p3d,r3d,v3d
#ifdef OCN_GISS_SM
     *    ,rx,ry,gx,gy,sx,sy
#endif
      USE OCEAN, only : OXYP,OLAT2D_DG,OJ_BUDG,OWTBUDG
      USE OCEAN, only : FOCEAN,HATMO,HOCEAN
      USE OCEAN, only : LMM,LMU,LMV
      USE OCEAN, only : nbyzmax,
     &     nbyzm,nbyzu,nbyzv,nbyzc,
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv, i1yzc,i2yzc

      USE OCEAN, only: alloc_odiff

      USE OCEAN_DYN, only : DH,BYDH,VBAR, dZGdP, GUP,GDN, SUP,SDN
      USE OCEAN_DYN, only : MMI,SMU,SMV,SMW,CONV,MU,MV,MW
      use Dictionary_mod, only : sync_param

#ifdef TRACERS_OCEAN
      use ocean, only : ntrtrans,asmu,asmv,asmw,motr,mosv0
#endif


      IMPLICIT NONE

      INTEGER :: IER
      integer :: img, jmg, lmg

C*
C**** Define the ocean (Russell) grid 
C*
      call init_oceanr_grid  
C****
 
      call getDomainBounds(ogrid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      call sync_param( "ocean_use_qus", use_qus )

      call alloc_straits

#ifdef TRACERS_OCEAN
      call alloc_ocn_tracer_com
#endif

      ALLOCATE(   LMM(IM,J_0H:J_1H), STAT = IER)
      ALLOCATE(   LMU(IM,J_0H:J_1H), STAT = IER)
      ALLOCATE(   LMV(IM,J_0H:J_1H), STAT = IER)
      lmm = 0
      lmu = 0
      lmv = 0

      ALLOCATE(   MO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   UO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   VO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   UOD(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   VOD(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( G0M (IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GXMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GYMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GZMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( S0M (IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SXMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SYMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SZMO(IM,J_0H:J_1H,LMO), STAT = IER)

#ifdef OCN_GISS_SM
      ALLOCATE( RX(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( RY(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GX(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GY(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SX(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SY(IM,J_0H:J_1H,LMO), STAT = IER)
#endif

      if(use_qus==1) then
      ALLOCATE( GXXMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GYYMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GZZMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GXYMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GYZMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GZXMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SXXMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SYYMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SZZMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SXYMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SYZMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SZXMO(IM,J_0H:J_1H,LMO), STAT = IER)
      gxxmo=0.; gyymo=0.; gzzmo=0.; gxymo=0.; gyzmo=0.; gzxmo=0.
      sxxmo=0.; syymo=0.; szzmo=0.; sxymo=0.; syzmo=0.; szxmo=0.
      endif

      if (am_i_root()) then
        img = im
        jmg = jm
        lmg = lmo
      else
        img = 1
        jmg = 1
        lmg = 1
      end if

      ALLOCATE( FOCEAN(IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( HATMO(IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( HOCEAN(IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OPRESS(IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OPBOT (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OGEOZ (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OGEOZ_SV (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( kpl   (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OXYP  (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OLAT2D_DG  (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OJ_BUDG  (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OWTBUDG  (IM,J_0H:J_1H), STAT = IER)

!!!   ALLOCATE(   PO(IM,J_0H:J_1H,LMO), STAT = IER)
!!!   ALLOCATE(  PHI(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   DH(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( BYDH(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( VBAR(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(dZGdP(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  GUP(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  GDN(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  SUP(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  SDN(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  MMI(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  SMU(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  SMV(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  SMW(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( CONV(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   MU(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   MV(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   MW(IM,J_0H:J_1H,LMO), STAT = IER)

C**** Necessary initiallisation?
      MU=0. ; MV=0. ; MW=0. ; CONV=0. ; MMI=0.
      UO=0. ; VO=0.
      SMU=0.; SMV=0.; SMW=0.; kpl=3

      allocate( g3d (lmo,im,j_0h:j_1h) )
      allocate( t3d (lmo,im,j_0h:j_1h) )
      allocate( s3d (lmo,im,j_0h:j_1h) )
      allocate( p3d (lmo,im,j_0h:j_1h) )
      allocate( r3d (lmo,im,j_0h:j_1h) )
      allocate( v3d (lmo,im,j_0h:j_1h) )
      g3d = 0.; t3d = 0.; s3d = 0.; p3d = 0.; r3d = 0.; v3d = 0.

      ALLOCATE(NBYZM(J_0H:J_1H,LMO))
      ALLOCATE(NBYZU(J_0H:J_1H,LMO))
      ALLOCATE(NBYZV(J_0H:J_1H,LMO))
      ALLOCATE(NBYZC(J_0H:J_1H,LMO))
      ALLOCATE(I1YZM(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I2YZM(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I1YZU(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I2YZU(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I1YZV(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I2YZV(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I1YZC(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I2YZC(NBYZMAX,J_0H:J_1H,LMO))

#ifdef TRACERS_OCEAN
      allocate( motr(im,j_0h:j_1h,lmo), stat = ier)
      motr = 0.
      call sync_param('ocean_ntrtrans',ntrtrans)
      if(ntrtrans.gt.1) then
#ifdef TRACERS_SPECIAL_O18
! Salt/freshwater are not (yet) transported on the ntrtrans schedule.  Full
! instantaneous consistency of water isotopes with freshwater requires
! identical timestepping.   Future releases will relax this restriction.
        call stop_model('ocean_ntrtrans>1 incompatible with '//
     &       'TRACERS_SPECIAL_O18',255)
#endif
        allocate( mosv0(im,j_0h:j_1h,lmo), stat = ier)
        allocate( asmu(im,j_0h:j_1h,lmo), stat = ier)
        allocate( asmv(im,j_0h:j_1h,lmo), stat = ier)
        allocate( asmw(im,j_0h:j_1h,lmo), stat = ier)
        mosv0 = 0.; asmu = 0.; asmv = 0.; asmw = 0.
      endif
#endif

c??   call ALLOC_GM_COM(agrid)
c      call ALLOC_KPP_COM(ogrid) ! alloc deferred until lsrpd known
#ifdef OCN_GISS_TURB
      call alloc_gissmix_com(ogrid)
#endif
#ifdef OCN_GISS_SM
      call alloc_giss_sm_com(ogrid)
#endif
      call alloc_odiag(ogrid)
      !call alloc_afluxes
      !call ALLOC_OFLUXES(atmocn)

#ifdef TRACERS_OceanBiology
      call alloc_obio_com
#endif
      call alloc_odiff(ogrid)

      call alloc_ocnmeso_com

      call read_ocean_topo
      if(ogrid%have_domain) CALL GEOMO

      return
      end subroutine alloc_ocean

      subroutine read_ocean_topo
C**** READ IN LANDMASKS AND TOPOGRAPHIC DATA
      USE OCEAN, only : FOCEAN,HATMO,HOCEAN
      use pario, only : par_open,par_close,read_dist_data
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : halo_update
      IMPLICIT NONE
      INTEGER :: fid

      fid = par_open(grid,'TOPO_OC','read')
      call read_dist_data(grid,fid,'focean',focean)
      call read_dist_data(grid,fid,'zatmo' ,hatmo )
      call read_dist_data(grid,fid,'zocean',hocean)
      call par_close(grid,fid)

      call halo_update(grid,focean)
      call halo_update(grid,hocean)
      call halo_update(grid,hatmo)

      return
      end subroutine read_ocean_topo
