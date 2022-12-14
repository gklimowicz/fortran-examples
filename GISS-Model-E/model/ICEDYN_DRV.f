c****
C**** ICEDYN_DRV.f    Sea ICE DYNamics    2006/12/21
C****
#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif

#ifdef CUBED_SPHERE
#ifdef HYCOM1degRefined
#define OCEAN_IMPORTEXPORT_ON_BGRID
#endif
#ifdef HYCOM1degUnrefined
#define OCEAN_IMPORTEXPORT_ON_BGRID
#endif
#endif

      MODULE ICEDYN_COM
!@sum  ICEDYN_COM holds global variables for dynamic sea ice
!@auth Gary Russell/Gavin Schmidt
      USE MDIAG_COM, only : lname_strlen,sname_strlen,units_strlen
#ifdef CUBED_SPHERE
      USE cs2ll_utils, only : cs2llint_type,ll2csint_type
#else
      use domain_decomp_1d, only : band_pack_type
#endif
      use cdl_mod
      USE EXCHANGE_TYPES, only : iceocn_xchng_vars
      IMPLICIT NONE
      SAVE

C**** Dimensions of ice advection grid (EDIT FOR ADVSI GRID CHANGE) are the same as
C**** atmospheric grid
      !INTEGER, parameter :: IMIC=IM, JMIC=JM

C**** Ice advection grid, same as atmospheric grid (CS or latlon)
C**** dimensions IMIC = IM, JMIC = JM
      !TYPE(DIST_GRID) :: grid_MIC

!@var igice
      type(iceocn_xchng_vars) :: igice

#ifdef CUBED_SPHERE
      type(cs2llint_type) :: CS2ICEint_a,CS2ICEint_b
      type(ll2csint_type) :: i2a_uc,i2a_vc ,ICE2CSint
c arrays for sea ice advection
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FOA,BYFOA,CONNECT
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: UVLLATUC,UVLLATVC
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: UVMATUC,UVMATVC
#else
! To allow the ice dynamics to run on a subset of processors:
!@var pack_a2i,pack_i2a contain info for redistributing data from
!@+   atmos. domains to icedyn domains and vice versa.
      type(band_pack_type), target :: pack_a2i,pack_i2a
#ifdef EXPEL_COASTAL_ICEXS
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: CONNECT
#endif
#endif

!@var aUSI,aVSI sea ice velocities (B-grid) having atmos. domain bounds,
!@+   needed in the lat-lon version of ADVSI
      !REAL*8, ALLOCATABLE, DIMENSION(:,:) :: aUSI,aVSI

C**** Ice dynamics diagnostics
      INTEGER, PARAMETER :: KICIJ=7
!@var IJ_xxx Names for ICIJ diagnostics
      INTEGER IJ_USI,IJ_VSI,IJ_DMUI,IJ_DMVI,IJ_PICE,IJ_RSI,IJ_SISPEED
!@var ICIJ lat-lon ice dynamic diagnostics (on ice dyn. grid)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: ICIJ!,ICIJg
!@var lname_icij Long names for ICIJ diagnostics
      CHARACTER(len=lname_strlen), DIMENSION(KICIJ) :: LNAME_ICIJ
!@var sname_icij Short names for ICIJ diagnostics
      CHARACTER(len=sname_strlen), DIMENSION(KICIJ) :: SNAME_ICIJ
!@var units_icij Units for ICIJ diagnostics
      CHARACTER(len=units_strlen), DIMENSION(KICIJ) :: UNITS_ICIJ
!@var ia_icij IDACC numbers for ICIJ diagnostics
      INTEGER, DIMENSION(KICIJ) :: IA_ICIJ
!@var scale_icij scales for ICIJ diagnostics
      REAL*8, DIMENSION(KICIJ) :: SCALE_ICIJ
!@var [ij]grid_icij Grid descriptors for ICIJ diagnostics
       INTEGER, DIMENSION(KICIJ) :: IGRID_ICIJ,JGRID_ICIJ
!@var denom_icij denominators for ICIJ diagnostics
       INTEGER, DIMENSION(KICIJ) :: DENOM_ICIJ

!@var cdl_icij consolidated metadata for ICIJ output fields in CDL notation
       type(cdl_type) :: cdl_icij

      END MODULE ICEDYN_COM

      SUBROUTINE ALLOC_ICEDYN_COM(grid_atm)
!@sum ALLOC_ICEDYN_COM allocates arrays defined in the ICEDYN_COM module
!@auth Rosalinda de Fainchtein

      USE DOMAIN_DECOMP_1D, only : getDomainBounds,DIST_GRID,am_I_root
      !USE ICEDYN_COM, only : grid_MIC
      USE ICEDYN_COM, only : KICIJ
      USE ICEDYN_COM, only : ICIJ,igice
#ifdef CUBED_SPHERE
      USE ICEDYN_COM, only : FOA,BYFOA,CONNECT
#else
#ifdef EXPEL_COASTAL_ICEXS
      USE ICEDYN_COM, only : CONNECT
#endif
#endif
      USE ICEDYN, only : grid_icdyn
      IMPLICIT NONE

      LOGICAL, SAVE :: init=.false.
      INTEGER :: I_1H    , I_0H
      INTEGER :: J_1H    , J_0H
      INTEGER :: I_1H_MIC, I_0H_MIC
      INTEGER :: J_1H_MIC, J_0H_MIC
      INTEGER :: IER
      TYPE(DIST_GRID) :: grid_atm

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

C**** Allocate arrays defined on the ice rheology grid

      call getDomainBounds(grid_icdyn,
     &     I_STRT_HALO=I_0H, I_STOP_HALO=I_1H,
     &     J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      igice%grid => grid_icdyn

      ALLOCATE(igice % DMUI( I_0H:I_1H , J_0H:J_1H ),
     &         igice % DMVI( I_0H:I_1H , J_0H:J_1H ))

      ALLOCATE(igice % UOSURF( I_0H:I_1H , J_0H:J_1H ),
     &         igice % VOSURF( I_0H:I_1H , J_0H:J_1H ))
      igice % UOSURF = 0. ! in case there is no dynamic ocean
      igice % VOSURF = 0. ! ""

      ALLOCATE(  ICIJ(I_0H:I_1H, J_0H:J_1H,KICIJ),
     &     STAT = IER)

c      if(am_I_root()) then
c         allocate(ICIJg(imic,jmic,KICIJ))
c      else
c         allocate(ICIJg(1,1,1))
c      end if


#ifdef CUBED_SPHERE
C**** Allocate ice advection arrays defined on the atmospheric grid
      !grid_MIC=grid_atm

      call getDomainBounds(grid_atm,
     &     I_STRT_HALO=I_0H_MIC, I_STOP_HALO=I_1H_MIC,
     &     J_STRT_HALO=J_0H_MIC, J_STOP_HALO=J_1H_MIC)

      ALLOCATE(   FOA(I_0H_MIC:I_1H_MIC, J_0H_MIC:J_1H_MIC),
     &          BYFOA(I_0H_MIC:I_1H_MIC, J_0H_MIC:J_1H_MIC),
     &        CONNECT(I_0H_MIC:I_1H_MIC, J_0H_MIC:J_1H_MIC))
#else
#ifdef EXPEL_COASTAL_ICEXS
      ALLOCATE(CONNECT(grid_atm%i_strt:grid_atm%i_stop,
     &                 grid_atm%j_strt_halo:grid_atm%j_stop_halo))
#endif
#endif

      return
      END SUBROUTINE ALLOC_ICEDYN_COM

      SUBROUTINE gather_icdiags ()
!@sum  collect the local acc-arrays into global arrays
!@+    run-time
!@auth Reto Ruedy

       USE ICEDYN_COM

       use domain_decomp_1d, only : pack_data

       IMPLICIT NONE

       END SUBROUTINE gather_icdiags

      subroutine def_rsf_icedyn(fid)
!@sum  def_rsf_icedyn defines ice dynam array structure in restart files
!@auth M. Kelley
!@ver  beta
      use icedyn, only : grid_icdyn,usi,vsi
      use icedyn_com, only : igice
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid_icdyn,fid,usi,'usi(dist_imic,dist_jmic)')
      call defvar(grid_icdyn,fid,vsi,'vsi(dist_imic,dist_jmic)')
      call defvar(grid_icdyn,fid,igice%uosurf,
     &     'uosurf_icdyn(dist_imic,dist_jmic)')
      call defvar(grid_icdyn,fid,igice%vosurf,
     &     'vosurf_icdyn(dist_imic,dist_jmic)')
      return
      end subroutine def_rsf_icedyn

      subroutine new_io_icedyn(fid,iaction)
!@sum  new_io_icedyn read/write ice dynam arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use icedyn, only : grid_icdyn,usi,vsi
      use icedyn_com, only : igice
      use model_com, only : ioread,iowrite
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to restart file
        call write_dist_data(grid_icdyn, fid, 'usi', usi)
        call write_dist_data(grid_icdyn, fid, 'vsi', vsi)
        call write_dist_data(grid_icdyn, fid, 'uosurf_icdyn',
     &       igice%uosurf)
        call write_dist_data(grid_icdyn, fid, 'vosurf_icdyn',
     &       igice%vosurf)
      case (ioread)            ! input from restart file
        call read_dist_data(grid_icdyn, fid, 'usi', usi)
        call read_dist_data(grid_icdyn, fid, 'vsi', vsi)
        call read_dist_data(grid_icdyn, fid, 'uosurf_icdyn',
     &       igice%uosurf)
        call read_dist_data(grid_icdyn, fid, 'vosurf_icdyn',
     &       igice%vosurf)
      end select
      return
      end subroutine new_io_icedyn

      subroutine def_rsf_icdiag(fid,r4_on_disk)
!@sum  def_rsf_icdiag defines ice diag array structure in restart/acc files
!@auth M. Kelley
!@ver  beta
      use icedyn, only : grid=>grid_icdyn
      use icedyn_com, only : icij
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
      call defvar(grid,fid,icij,'icij(dist_imic,dist_jmic,kicij)',
     &     r4_on_disk=r4_on_disk)
      return
      end subroutine def_rsf_icdiag

      subroutine new_io_icdiag(fid,iaction)
!@sum  new_io_icdiag read/write ice diag arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use icedyn, only : grid=>grid_icdyn
      use icedyn_com, only : icij
      use model_com, only : ioread,iowrite
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to restart or acc file
        call write_dist_data(grid, fid, 'icij', icij)
      case (ioread)            ! input from restart or acc file
        call read_dist_data(grid, fid, 'icij', icij)
      end select
      return
      end subroutine new_io_icdiag

      subroutine def_meta_icdiag(fid)
!@sum  def_meta_icdiag defines icedyn metadata in acc files
!@auth M. Kelley
!@ver  beta
      use icedyn, only : grid=>grid_icdyn
      use icedyn_com, only :
     &     ia_icij,denom_icij,scale_icij,sname_icij,cdl_icij
      use pario, only : defvar,write_attr
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_attr(grid,fid,'icij','reduction','sum')
      call write_attr(grid,fid,'icij','split_dim',3)
      call defvar(grid,fid,ia_icij,'ia_icij(kicij)')
      call defvar(grid,fid,denom_icij,'denom_icij(kicij)')
      call defvar(grid,fid,scale_icij,'scale_icij(kicij)')
      call defvar(grid,fid,sname_icij,'sname_icij(sname_strlen,kicij)')
      call defvar_cdl(grid,fid,cdl_icij,
     &     'cdl_icij(cdl_strlen,kcdl_icij)')

      return
      end subroutine def_meta_icdiag

      subroutine write_meta_icdiag(fid)
!@sum  write_meta_icdiag write icedyn accumulation metadata to file
!@auth M. Kelley
!@ver  beta
      use icedyn, only : grid=>grid_icdyn
      use icedyn_com, only :
     &     ia_icij,denom_icij,scale_icij,sname_icij,cdl_icij
      use pario, only : write_dist_data,write_data
      use cdl_mod, only : write_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_data(grid,fid,'ia_icij',ia_icij)
      call write_data(grid,fid,'denom_icij',denom_icij)
      call write_data(grid,fid,'scale_icij',scale_icij)
      call write_data(grid,fid,'sname_icij',sname_icij)
      call write_cdl(grid,fid,'cdl_icij',cdl_icij)

      return
      end subroutine write_meta_icdiag

      subroutine set_ioptrs_iceacc_default
c point i/o pointers for diagnostic accumlations to the
c instances of the arrays used during normal operation.
c temporarily empty.
      return
      end subroutine set_ioptrs_iceacc_default

      SUBROUTINE reset_icdiag
!@sum reset_icdiag resets ice dynamic diagnostic arrays
!@auth Gavin Schmidt
      USE ICEDYN_COM
      IMPLICIT NONE

      ICIJ=0.
      RETURN
      END SUBROUTINE reset_icdiag

      SUBROUTINE DYNSI(atmice,iceocn,si_ocn)
!@sum  DYNSI calculates ice velocites
!@+    Note that the ice velocities are calculated on the ice grid
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang), D. Gueyffier (cubed sphere)
!@auth M. Kelley (cubed sphere)
      USE CONSTANT, only : rhoi,grav,omega,rhows
      USE MODEL_COM, only : dts=>dtsrc
      USE ICEDYN, only : imicdyn,jmicdyn,  !dimensions of icedyn grid
     &     nx1,ny1
      USE DOMAIN_DECOMP_1D, only : ICE_HALO=>HALO_UPDATE,
     &     NORTH, SOUTH, hasSouthPole, hasNorthPole, getDomainBounds
      USE ICEDYN, only : dxyn,dxys,bydxyp,dxyv,dxp,dyv
      USE ICEDYN, only : grid_ICDYN
      USE ICEDYN, only : press,heffm,uvm,dwatn,cor
     *     ,sinwat,coswat,bydts,sinen,uice,vice,heff,area,gairx,gairy
     *     ,gwatx,gwaty,pgfub,pgfvb,amass,dmu,dmv
     *     ,usi,vsi,iFOCEAN=>FOCEAN
      USE ICEDYN_COM, only : icij,ij_usi
     *     ,ij_vsi,ij_dmui,ij_dmvi,ij_pice,ij_rsi,igice,ij_sispeed
      USE SEAICE, only : ace1i
#ifdef CUBED_SPHERE
      use icedyn_com, only : CS2ICEint_a,CS2ICEint_b,ICE2CSint
      use cs2ll_utils, only : cs2llint_lij,cs2llint_lluv
      use cs2ll_utils, only : ll2csint_ij
#else
      USE DOMAIN_DECOMP_1D, only : band_pack,
     &     hasSouthPole, hasNorthPole
      !needs evaluation USE GEOM, only : cosip,sinip
      use icedyn_com, only : pack_a2i,pack_i2a
#endif
      USE TimerPackage_mod, only: startTimer => start
      USE TimerPackage_mod, only: stopTimer => stop
      USE EXCHANGE_TYPES, only : atmice_xchng_vars,iceocn_xchng_vars
      USE SEAICE_COM, only : icestate
      IMPLICIT NONE
      type(atmice_xchng_vars) :: atmice
      type(iceocn_xchng_vars) :: iceocn
      type(icestate) :: si_ocn
      SAVE
C**** intermediate calculation for pressure gradient terms
      REAL*8, allocatable, DIMENSION(:,:) :: pgfu, pgfv
C****
      real*8, allocatable, dimension(:,:) ::
     &     aPtmp,iPtmp,iRSI,iMSI,iDMUA,iDMVA,iUI2rho,admu,admv
      real*8, allocatable, dimension(:,:,:) :: alij_tmp,ilij_tmp

      REAL*8, DIMENSION(:,:), POINTER ::
     &     RSI,MSI,SNOWI,FOCEAN,RSISAVE,DMUA,DMVA,
     &     UOSURF,VOSURF,DMUI,DMVI

      REAL*8, PARAMETER :: BYRHOI=1D0/RHOI
      REAL*8 :: hemi
      INTEGER :: I,J,ip1,im1
      REAL*8 :: DMUINP,duA,dvA,rsib,USISPEED,VSISPEED
      INTEGER :: aIM, aJM
      INTEGER :: iJ_1   , iJ_0
      INTEGER :: iJ_1S  , iJ_0S
      INTEGER :: iJ_1H  , iJ_0H
      INTEGER :: iJ_1STG, iJ_0STG
      INTEGER :: aI_1   , aI_0
      INTEGER :: aJ_1   , aJ_0
      INTEGER :: aI_1H  , aI_0H
      INTEGER :: aJ_1H  , aJ_0H
      INTEGER :: aJ_1S  , aJ_0S

      call startTimer('DYNSI()')

      FOCEAN => ATMICE%FOCEAN
      DMUA => ATMICE%DMUA
      DMVA => ATMICE%DMVA
      RSI => SI_OCN%RSI
      MSI => SI_OCN%MSI
      SNOWI => SI_OCN%SNOWI
      RSISAVE => SI_OCN%RSISAVE
      UOSURF => IGICE%UOSURF
      VOSURF => IGICE%VOSURF
      DMUI => IGICE%DMUI
      DMVI => IGICE%DMVI

C**** Get loop indices  corresponding to grid_ICDYN and atm. grid structures
      call getDomainBounds(grid_ICDYN, J_STRT=iJ_0, J_STOP=iJ_1,
     &     J_STRT_SKP=iJ_0S   , J_STOP_SKP=iJ_1S,
     &     J_STRT_HALO=iJ_0H  , J_STOP_HALO=iJ_1H )
      call getDomainBounds(grid_ICDYN, J_STRT_STGR=iJ_0STG,
     &     J_STOP_STGR=iJ_1STG)
      call getDomainBounds(atmice%grid, I_STRT=aI_0, I_STOP=aI_1,
     &     J_STRT=aJ_0, J_STOP=aJ_1)
      call getDomainBounds(atmice%grid,
     &     I_STRT_HALO=aI_0H, I_STOP_HALO=aI_1H,
     &     J_STRT_HALO=aJ_0H  , J_STOP_HALO=aJ_1H )
      call getDomainBounds(atmice%grid,
     &     J_STRT_SKP=aJ_0S, J_STOP_SKP=aJ_1S)
      aIM = atmice%grid%im_world
      aJM = atmice%grid%jm_world

      allocate(
     &     pgfu(IMICDYN,
     &     grid_ICDYN%J_STRT_HALO:grid_ICDYN%J_STOP_HALO),
     &     pgfv(IMICDYN,
     &     grid_ICDYN%J_STRT_HALO:grid_ICDYN%J_STOP_HALO)
     &     )
      allocate(
     &     aPtmp(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     iPtmp(1:IMICDYN,iJ_0H:iJ_1H),
     &     iRSI(1:IMICDYN,iJ_0H:iJ_1H),
     &     iMSI(1:IMICDYN,iJ_0H:iJ_1H)
     &     )

C**** Start main loop
C**** Replicate polar boxes
      if (hasNorthPole(atmice%grid)) then
        RSI(2:aIM,aJM)=RSI(1,aJM)
        MSI(2:aIM,aJM)=MSI(1,aJM)
        SNOWI(2:aIM,aJM)=SNOWI(1,aJM)
        DMUA(2:aIM,aJM) = DMUA(1,aJM)
        DMVA(2:aIM,aJM) = DMVA(1,aJM)
      end if

c**** interpolate air stress from A grid in atmos, to B grid in ice
C**** change of unit from change of momentum, to flux

      DO J=aJ_0,aJ_1
        do i=aI_0,aI_1
          IF (FOCEAN(I,J)*RSI(I,J).le.0) THEN
            DMUA(I,J) = 0.
            DMVA(I,J) = 0.
          END IF
        END DO
      END DO

c**** getting instance of (DMUA, DVMA) on the icedyn grid
      allocate(iDMUA(1:IMICDYN,iJ_0H:iJ_1H),
     &         iDMVA(1:IMICDYN,iJ_0H:iJ_1H))
#ifdef CUBED_SPHERE
      call cs2llint_lluv(atmice%grid,CS2ICEint_b,dmua,dmva,
     &     idmua,idmva)
      do j=iJ_0,iJ_1S
        do i=1,imicdyn
          gairx(i+1,j) = idmua(i,j)*bydts
          gairy(i+1,j) = idmva(i,j)*bydts
        enddo
        gairx((/1,nx1/),j) = gairx((/nx1-1,2/),j)
        gairy((/1,nx1/),j) = gairy((/nx1-1,2/),j)
      enddo
#else
      call band_pack(pack_a2i, DMUA, iDMUA) ! fills halo
      call band_pack(pack_a2i, DMVA, iDMVA) ! fills halo

c needs evaluation      if (grid_icdyn%HAVE_NORTH_POLE) then
c needs evaluation        dua = idmua(1,jmicdyn)
c needs evaluation        dva = idmva(1,jmicdyn)
c needs evaluation        idmua(:,jmicdyn)=dua*cosip(:)+dva*sinip(:)
c needs evaluation        idmva(:,jmicdyn)=dva*cosip(:)-dua*sinip(:)
c needs evaluation      end if
      do j=iJ_0,iJ_1S
        im1=imicdyn
        do i=1,imicdyn
          GAIRX(i,j)=0.25*(idmua(i,j)+idmua(im1,j)
     &                    +idmua(im1,j+1)+idmua(i,j+1))*bydts
          GAIRY(i,j)=0.25*(idmva(i,j)+idmva(im1,j)
     &                    +idmva(im1,j+1)+idmva(i,j+1))*bydts
          im1=i
        enddo
      enddo
      IF (hasNorthPole(grid_ICDYN)) THEN
        GAIRX(1:nx1,jmicdyn)=idmua(1,jmicdyn)*bydts
        GAIRY(1:nx1,jmicdyn)=idmva(1,jmicdyn)*bydts
      END IF
      do j=iJ_0,iJ_1S
       GAIRX(nx1-1,j)=GAIRX(1,j)
       GAIRY(nx1-1,j)=GAIRY(1,j)
       GAIRX(nx1,j)=GAIRX(2,j)
       GAIRY(nx1,j)=GAIRY(2,j)
      enddo
#endif
      deallocate(iDMUA,iDMVA)

C**** save current value of sea ice concentration for ADVSI
C**** RSISAVE is on atmospheric grid
      RSISAVE(:,:)=RSI(:,:)

C**** Pressure anomaly at surface APRESS: calculated by sea ice routines
C**** APRESS is on atmospheric grid. We are now no longer using this as
C**** a forcing in the sea ice dynamics because it is partially
C**** included already in the internal ice pressure gradient. The
C**** atmospheric pressure gradient term does not in general produce a
C**** horizontal force in a solid (such as ice).

C**** calculate sea surface tilt on atmospheric C grid
C**** (using OGEOZA on atmospheric grid plus displacement of free
C**** surface due to presence of ice). This is ignored in favour of
C**** geostrophy if osurf_tilt=0.
C**** PGF is an accelaration

C****  define scalar pressure on atm grid then regrid it to the icedyn grid
      DO J=aJ_0,aJ_1
         DO I=aI_0,aI_1
             aPtmp(I,J)=(iceocn%OGEOZA(I,J)
     *            +(RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I))*GRAV/RHOWS)
        END DO
      END DO

#ifdef CUBED_SPHERE
c bundle the qtys to be interpolated
      allocate(alij_tmp(3,aI_0H:aI_1H,aJ_0H:aJ_1H))
      allocate(ilij_tmp(3,1:IMICDYN,iJ_0H:iJ_1H))
      do j=aJ_0,aJ_1
        do i=aI_0,aI_1
          alij_tmp(1,i,j) = aPtmp(i,j)
          alij_tmp(2,i,j) = RSI(i,j)
          alij_tmp(3,i,j) = MSI(i,j)
        enddo
      enddo
      call cs2llint_lij(atmice%grid,CS2ICEint_a,alij_tmp,ilij_tmp)
      do j=iJ_0,iJ_1
        do i=1,IMICDYN
          iPtmp(i,j) = ilij_tmp(1,i,j)
          iRSI(i,j)  = ilij_tmp(2,i,j)
          iMSI(i,j)  = ilij_tmp(3,i,j)
        enddo
      enddo
      deallocate(alij_tmp,ilij_tmp)
#else
      call band_pack(pack_a2i, aPtmp, iPtmp) ! like iPtmp = aPtmp
      call band_pack(pack_a2i, RSI, iRSI)    !      iRSI = RSI
      call band_pack(pack_a2i, MSI, iMSI)    !      iMSI = MSI
#endif


c-------------------------------------------------------------------
c Begin icedyn-processors-only code region
      icedyn_processors_only: if(grid_ICDYN%have_domain) then
c-------------------------------------------------------------------

      CALL ICE_HALO(grid_ICDYN, iPtmp , from=NORTH )
      CALL ICE_HALO(grid_ICDYN, iRSI   , from=NORTH )

c*** Calculate gradient on ice dyn. grid
      if (hasNorthPole(grid_ICDYN)) PGFU(1:IMICDYN,JMICDYN)=0
      if (hasSouthPole(grid_ICDYN)) PGFU(1:IMICDYN, 1)=0  !RKF
      DO J=iJ_0S,iJ_1S
        I=IMICDYN
        DO IP1=1,IMICDYN
          IF(iFOCEAN(I,J).gt.0 .and. iFOCEAN(IP1,J).gt.0. .and.
     *         iRSI(I,J)+iRSI(IP1,J).gt.0.) THEN
            PGFU(I,J)=-(iPtmp(IP1,J)-iPtmp(I,J))/DXP(J)
          ELSE
            PGFU(I,J)=0.
          END IF
          I=IP1
        END DO
      END DO

      DO J=iJ_0,iJ_1S
        DO I=1,IMICDYN
          IF(iFOCEAN(I,J+1).gt.0 .and. iFOCEAN(I,J).gt.0. .and.
     *         iRSI(I,J)+iRSI(I,J+1).gt.0.) THEN
            PGFV(I,J)=-(iPtmp(I,J+1)-iPtmp(I,J))/DYV(J+1)
          ELSE
            PGFV(I,J)=0.
          END IF
        END DO
      END DO
c*

C**** Set up ice grid variables
C**** HEFF,AREA on primary (tracer) grid for ice
      do j=iJ_0,iJ_1
         do i=2,NX1-1
            HEFF(I,J)=iRSI(I-1,J)*(ACE1I+iMSI(I-1,J))*BYRHOI
            HEFF(I,J)=HEFF(I,J)*HEFFM(I,J)
            AREA(I,J)=iRSI(I-1,J)
        enddo
      enddo

C**** fill in overlap regions
      DO J=iJ_0,iJ_1
        HEFF(1,J)=HEFF(NX1-1,J)
        AREA(1,J)=AREA(NX1-1,J)
        HEFF(NX1,J)=HEFF(2,J)
        AREA(NX1,J)=AREA(2,J)
      END DO

C****
C**** Set up mass per unit area and coriolis term (on ice grid B)
C****
C**** Update halo for HEFF
      CALL ICE_HALO(grid_ICDYN, HEFF, from=NORTH    )

      DO J=iJ_0,iJ_1S
      DO I=1,NX1-1
        AMASS(I,J)=RHOI*0.25*(HEFF(I,J)
     *       +HEFF(I+1,J)+HEFF(I,J+1)+HEFF(I+1,J+1))
        COR(I,J)=AMASS(I,J)*2.0*OMEGA*SINEN(I,J)
      END DO
      END DO
c**** set north pole
      if (hasNorthPole(grid_ICDYN)) then
        do i=1,nx1
          AMASS(i,jmicdyn)= RHOI*HEFF(1,JMICDYN)
          COR  (i,jmicdyn)= AMASS(i,jmicdyn)
     *         *2.0*OMEGA*SINEN(1,JMICDYN)
        end do
      end if                    !end NORTH_POLE block if

c**** interpolate air, current and ice velocity from C grid to B grid
C**** This should be more generally from ocean grid to ice grid
C**** NOTE: UOSURF, VOSURF are expected to be on the A-grid

C**** Update halo for USI,UOSURF,VOSURF,PGFU
      CALL ICE_HALO(grid_ICDYN, USI   , from=NORTH    )
      CALL ICE_HALO(grid_ICDYN, UOSURF , from=NORTH    )
      CALL ICE_HALO(grid_ICDYN, VOSURF , from=NORTH    )
      CALL ICE_HALO(grid_ICDYN, PGFU  , from=NORTH    )

      do j=iJ_0,iJ_1S
        im1=imicdyn
        do i=1,imicdyn
#ifdef OCEAN_IMPORTEXPORT_ON_BGRID
          GWATX(i,j)=UOSURF(im1,j)
          GWATY(i,j)=VOSURF(im1,j)
#else
          GWATX(i,j)=0.25*(UOSURF(im1,j)  +UOSURF(im1,j+1)
     &         +UOSURF(i,j)+UOSURF(i,j+1))                     ! ocean -> iceB
          GWATY(i,j)=0.25*(VOSURF(im1,j)  +VOSURF(im1,j+1)
     &         +VOSURF(i,j)+VOSURF(i,j+1))                     ! y component
#endif
          PGFUB(i,j)=0.5*(PGFU(im1,j)  +PGFU(im1,j+1))   ! iceC--> iceB
          PGFVB(i,j)=0.5*(PGFV(im1,j)  +PGFV(i,j))       ! y component
          im1=i
        enddo
      enddo
c**** set north pole
      if (hasNorthPole(grid_ICDYN)) then
        do i=1,imicdyn
          GWATX(i,jmicdyn)=0.
          PGFUB(i,jmicdyn)=0.
          GWATY(i,jmicdyn)=0.
          PGFVB(i,jmicdyn)=0.
        enddo
      end if                    !end NORTH_POLE block if
      DO J=iJ_0,iJ_1
        GWATX(nx1-1,J)=GWATX(1,J)
        GWATY(nx1-1,J)=GWATY(1,J)
        PGFUB(nx1-1,J)=PGFUB(1,J)
        PGFVB(nx1-1,J)=PGFVB(1,J)
        GWATX(nx1,J)=GWATX(2,J)
        GWATY(nx1,J)=GWATY(2,J)
        PGFUB(nx1,J)=PGFUB(2,J)
        PGFVB(nx1,J)=PGFVB(2,J)
      END DO

c**** read in sea ice velocity
      DO J=iJ_0,iJ_1
        DO I=1,IMICDYN
          UICE(I+1,J,1)=USI(I,J)
          VICE(I+1,J,1)=VSI(I,J)
        END DO
        UICE(1  ,J,1)=USI(IMICDYN,J)
        UICE(NX1,J,1)=USI(1,      J)
        VICE(1  ,J,1)=VSI(IMICDYN,J)
        VICE(NX1,J,1)=VSI(1,      J)
        DO I=1,NX1
          UICE(I,J,2)=0.
          VICE(I,J,2)=0.
          UICE(I,J,3)=0.
          VICE(I,J,3)=0.
        END DO
      END DO

C**** do the looping over pseudo-timesteps
      CALL VPICEDYN

C**** Calculate stress on ice velocity grid (B grid)
      DO J=iJ_0,iJ_1
        hemi=-1.
        if (J.gt.NY1/2) hemi=1.
        DO I=1,NX1-1
          DMU(i,j)=DTS*dwatn(i,j)*(COSWAT*(UICE(i,j,1)-GWATX(i,j))-
     *         HEMI*SINWAT*(VICE(i,j,1)-GWATY(i,j)))
          DMV(i,j)=DTS*dwatn(i,j)*(HEMI*SINWAT*(UICE(i,j,1)-GWATX(i,j))
     *         +COSWAT*(VICE(i,j,1)-GWATY(i,j)))
        END DO
      END DO
      DO J=iJ_0,iJ_1
        DMU(1,J)=DMU(NX1-1,J)
        DMU(NX1,J)=DMU(2,J)
      END DO

C**** Save ice velocities in USI,VSI arrays
C**** Interpolate ice stress from its B grid to C grid
C**** Update halos for UICE and DMU
      CALL ICE_HALO(grid_ICDYN,  UICE, from=SOUTH     )
      CALL ICE_HALO(grid_ICDYN,   DMU, from=SOUTH     )

      do j=iJ_0S,iJ_1S
        do i=1,imicdyn
          usi(i,j)=uice(i+1,j,1)
          if(abs(usi(i,j)).lt.1d-10) usi(i,j)=0
          vsi(i,j)=vice(i+1,j,1)
          if(abs(vsi(i,j)).lt.1d-10) vsi(i,j)=0
        enddo
        i=imicdyn
        do ip1=1,imicdyn
C**** Rescale DMUI,DMVI to be net momentum into ocean
#ifdef OCEAN_IMPORTEXPORT_ON_BGRID
          rsib = (DXYN(J)*(
     &         iFOCEAN(I,J)*iRSI(I,J)+iFOCEAN(ip1,J)*iRSI(ip1,J)
     &         )  + DXYS(J+1)*(
     &       iFOCEAN(I,J+1)*iRSI(I,J+1)+iFOCEAN(ip1,J+1)*iRSI(ip1,J+1)
     &         ) )/DXYV(J+1)
          DMUI(I,J) = dmu(i+1,j)*rsib
          DMVI(I,J) = dmv(i+1,j)*rsib
#else
          DMUI(I,J) = 0.5*(dmu(i+1,j-1)+dmu(i+1,j))
          DMUI(I,J) = 0.5*DMUI(I,J)*
     &         (iFOCEAN(I,J)*iRSI(I,J)+iFOCEAN(ip1,J)*iRSI(ip1,J))
          DMVI(I,J) = 0.5*(dmv(i,j)+dmv(i+1,j))
          DMVI(I,J) = 0.5*DMVI(I,J)*
     &          (iFOCEAN(I,J  )*iRSI(I,J  )*DXYN(J)
     &          +iFOCEAN(I,J+1)*iRSI(I,J+1)*DXYS(J+1))
     &          /DXYV(J+1)
#endif
          i=ip1
        enddo
      enddo

C**** set south pole
      if (hasSouthPole(grid_ICDYN)) then
        dmui(:,1)=0.
      endif

C**** set north pole
      IF (hasNorthPole(grid_ICDYN)) THEN
        USI(:,jmicdyn)=0.
        VSI(:,jmicdyn)=0.
#ifndef OCEAN_IMPORTEXPORT_ON_BGRID
        j=jmicdyn
        i=imicdyn
        do ip1=1,imicdyn
          DMUI(I,J) = 0.5*(dmu(i+1,j-1)+dmu(i+1,j))
          DMUI(I,J) = 0.5*DMUI(I,J)*
     &         (iFOCEAN(I,J)*iRSI(I,J)+iFOCEAN(ip1,J)*iRSI(ip1,J))
          i=ip1
        enddo
#endif
        DMUINP=0.
        do i=1,imicdyn
          DMUINP = DMUINP + DMUI(i,jmicdyn)
        enddo
        DMUINP=DMUINP/imicdyn
        DMUI(:,jmicdyn)=DMUINP
        DMVI(:,jmicdyn)=0.
      END IF

c*** diagnostics
      DO J=iJ_0,iJ_1S
        DO I=1,imicdyn
         ip1=i+1
         USISPEED=0.; VSISPEED=0.
         if (ip1 .eq. IMICDYN+1) ip1=1
          IF (iFOCEAN(I,J).gt.0 .and. iFOCEAN(IP1,J).gt.0. .and.
     *         iRSI(I,J)+iRSI(IP1,J).gt.1d-4) THEN
            ICIJ(I,J,IJ_USI) =ICIJ(I,J,IJ_USI) +(iRSI(I,J)+iRSI(IP1,J))
     *           *0.5*(uice(i+1,j-1,1)+uice(i+1,j,1))
            ICIJ(I,J,IJ_DMUI)=ICIJ(I,J,IJ_DMUI)+DMUI(i,j)
            USISPEED=(iRSI(I,J)+iRSI(IP1,J))*0.5
     *              *(uice(i+1,j-1,1)+uice(i+1,j,1))
          END IF
          IF (iFOCEAN(I,J+1).gt.0 .and. iFOCEAN(I,J).gt.0. .and.
     *         iRSI(I,J)+iRSI(I,J+1).gt.1d-4) THEN
            ICIJ(I,J,IJ_VSI) =ICIJ(I,J,IJ_VSI) +(iRSI(I,J)+iRSI(I,J+1))
     *           *0.5*(vice(i,j,1)+vice(i+1,j,1))
            ICIJ(I,J,IJ_DMVI)=ICIJ(I,J,IJ_DMVI)+DMVI(i,j)
            VSISPEED=(iRSI(I,J)+iRSI(I,J+1))
     *           *0.5*(vice(i,j,1)+vice(i+1,j,1))
          END IF
C**** SISPEED: Speed of ice (i.e. mean absolute velocity)
C**** to account for back-and-forth movement of the ice
          ICIJ(I,J,IJ_SISPEED)=ICIJ(I,J,IJ_SISPEED)
     *      +sqrt(USISPEED**2+VSISPEED**2)
          ICIJ(I,J,IJ_PICE)=ICIJ(I,J,IJ_PICE)+ iRSI(I,J)*press(i+1,j)
          ICIJ(I,J,IJ_RSI) =ICIJ(I,J,IJ_RSI) + iRSI(I,J)
        END DO
      END DO
      IF (hasNorthPole(grid_ICDYN)) THEN
        ICIJ(1,JMICDYN,IJ_DMUI)=ICIJ(1,JMICDYN,IJ_DMUI)+DMUI(1,JMICDYN)
        ICIJ(1,JMICDYN,IJ_RSI) =ICIJ(1,JMICDYN,IJ_RSI) +iRSI(1,JMICDYN)
        ICIJ(1,JMICDYN,IJ_PICE)=ICIJ(1,JMICDYN,IJ_PICE)
     &       +iRSI(1,JMICDYN)*press(1,JMICDYN)
      END IF


c-------------------------------------------------------------------
c End icedyn-processors-only code region
      endif icedyn_processors_only
c-------------------------------------------------------------------

C**** Calculate ustar*2*rho for ice-ocean fluxes on atmosphere grid
C**** UI2rho = | tau |

#ifdef CUBED_SPHERE /* calculate stress magnitude before regrid */
      allocate(iUI2rho(1:IMICDYN,iJ_0H:iJ_1H))
      do j=iJ_0S,iJ_1S
        do i=1,imicdyn
          iUI2rho(i,j)= sqrt(dmu(i+1,j)**2 + dmv(i+1,j)**2) * bydts
        enddo
      enddo
      IF(hasNorthPole(grid_ICDYN)) iUI2rho(:,jmicdyn)=0.
      call ll2csint_ij(grid_icdyn,ICE2CSint,iUI2rho,iceocn%UI2rho)
      deallocate(iUI2rho)
      do j=aJ_0,aJ_1
        do i=aI_0,aI_1
          if(iceocn%FWATER(I,J)*RSI(i,j).eq.0) iceocn%UI2rho(i,j)=0
        enddo
      enddo
#else /* calculate stress magnitude after regrid */
C**** calculate 4 point average of B grid values of stresses
      allocate(admu(nx1,aJ_0H:aJ_1H))
      allocate(admv(nx1,aJ_0H:aJ_1H))
      call band_pack(pack_i2a, dmu, admu) ! fills halo
      call band_pack(pack_i2a, dmv, admv) ! fills halo
      do j=aJ_0,aJ_1
        do i=1,aim
          iceocn%UI2rho(i,j)=0
          if (iceocn%FWATER(I,J)*RSI(i,j).gt.0) THEN
            duA = 0.5*(DXYN(J)*(admu(i+1,j)+admu(i,j))
     *           +DXYS(j)*(admu(i+1,j-1)+admu(i,j-1)))*BYDXYP(J)
            dvA = 0.5*(DXYN(J)*(admv(i+1,j)+admv(i,j))
     *           +DXYS(j)*(admv(i+1,j-1)+admv(i,j-1)))*BYDXYP(J)
            iceocn%UI2rho(i,j)= sqrt (duA**2 + dvA**2) * bydts
          endif
        enddo
      enddo
      deallocate(admu,admv)
#endif

C**** Set uisurf,visurf (on atm A grid) for use in atmos. drag calc.
C**** uisurf/visurf are on atm grid but are latlon oriented
      call get_uisurf(usi,vsi,atmice)

      deallocate(pgfu, pgfv)
      deallocate(aPtmp,iPtmp,iRSI,iMSI)

      call stopTimer('DYNSI()')
      RETURN
      END SUBROUTINE DYNSI

#ifndef CUBED_SPHERE
      SUBROUTINE ADVSI(atmice)
!@sum  ADVSI advects sea ice
!@+    Currently set up to advect ice on AGCM grid (i.e. usidt/vsidt are
!@+    on the AGCM grid, and RSI/MSI/HSI etc. are unchanged)
!@+    At some point this will change (USIDT/VSIDT on ice grid, and RSI
!@+    etc. will need to be interpolated back and forth).
!@auth Gary Russell/Gavin Schmidt
      USE CONSTANT, only : grav,tf
      USE MODEL_COM, only :  kocean,dts=>dtsrc
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : SOUTH, NORTH
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE_COLUMN,
     &     hasNorthPole, hasSouthPole
      USE ICEDYN, only : dxyp,dyp,dxp,dxv,bydxyp
      !USE ICEDYN_COM, only : grid_MIC
      USE SEAICE, only : ace1i,xsi,Ti,
     * get_snow_ice_layer,set_snow_ice_layer,relayer,relayer_12
      USE SEAICE_COM, only : si_ocn,lmi
#ifdef TRACERS_WATER
      USE SEAICE, only : ntm
#endif
      USE EXCHANGE_TYPES, only : atmice_xchng_vars
#ifdef EXPEL_COASTAL_ICEXS
      USE ICEDYN_COM, only : connect
#endif
      IMPLICIT NONE
      type(atmice_xchng_vars) :: atmice

!@var USIDT,VSIDT sea ice fluxes, saved for advection (m)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: FMSI,
     &     FASI, FXSI, FYSI, FAW, BYFOA, USIDT, VSIDT
!@var NTRICE max. number of tracers to be advected (mass/heat/salt+)
      INTEGER :: NTRICE
#ifdef TRACERS_WATER
      INTEGER ITR
      REAL*8 TRSNOW(NTM,2), TRICE(NTM,LMI)
#endif
      REAL*8, ALLOCATABLE :: SFMSI(:),AMSI(:)
!@var MHS mass/heat/salt content of sea ice
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: FMSJ,MHS
      INTEGER I,J,L,IM1,IP1,K
      REAL*8 SFASI,DMHSI,ASI,YRSI,XRSI,FRSI,SICE,TMP,TICE,ENRG
      REAL*8 :: MICE(LMI),SICEg(LMI),HICE(LMI),TSIL(LMI),MSI1,FMSI2
      REAL*8 :: HSNOW(2),SNOWL(2),TSNW(2)

C****
C**** FLUXCB  USIDT  U compon of time integrated sea ice velocity (m)
C****         VSIDT  V compon of time integrated sea ice velocity (m)
C****
C****         FAW    flux of surface water area (m^2) = USIDT*DYP
C****         FASI   flux of sea ice area (m^2) = USIDT*DYP*RSIedge
C****         FMSI   flux of sea ice mass (kg) or heat (J) or salt (kg)

      REAL*8, DIMENSION(:,:), POINTER ::
     &     RSI,MSI,SNOWI,FOCEAN,AUSI,AVSI,RSIX,RSIY,RSISAVE
      REAL*8, DIMENSION(:,:,:), POINTER :: HSI,SSI
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:,:), POINTER :: TRSI
#endif

#ifdef EXPEL_COASTAL_ICEXS
!@var coastfac[xy]: A proportionality factor to compute the component
!@+   of advective velocity which limits ice buildup along
!@+   coastlines.  (At some gridcells, negative feedbacks on
!@+   ice production are not able to assert themselves when
!@+   sea ice does not reside on the ocean grid.)
      real*8 :: coastfacx,coastfacy
      real*8 :: du,dv
      integer :: cxi,cxip1,cyj,cyjp1
#endif

      INTEGER :: IM,JM
      INTEGER J_0, J_1, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
C**** Get grid parameters
      IM = si_ocn%grid%IM_WORLD
      JM = si_ocn%grid%JM_WORLD
      call getDomainBounds(si_ocn%grid, J_STRT=J_0,    J_STOP=J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S ,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H ,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)


#ifndef TRACERS_WATER
      NTRICE=3*(LMI+2)
#else
      NTRICE=(3+NTM)*(LMI+2)
#endif

      ALLOCATE(SFMSI(NTRICE),AMSI(NTRICE))
      ALLOCATE(FMSI(NTRICE,IM),FMSJ(IM,NTRICE,J_0H:J_1H))
      ALLOCATE(MHS(NTRICE,IM,J_0H:J_1H))
      ALLOCATE(FASI(IM,J_0H:J_1H)
     &     ,FXSI(IM,J_0H:J_1H)
     &     ,FYSI(IM,J_0H:J_1H)
     &     ,FAW(IM,J_0H:J_1H)
     &     ,BYFOA(IM,J_0H:J_1H) )

      ALLOCATE( USIDT(IM, J_0H:J_1H), VSIDT(IM, J_0H:J_1H) )


      FOCEAN => ATMICE%FOCEAN
      RSI => SI_OCN%RSI
      RSIX => SI_OCN%RSIX
      RSIY => SI_OCN%RSIY
      RSISAVE => SI_OCN%RSISAVE
      MSI => SI_OCN%MSI
      HSI => SI_OCN%HSI
      SSI => SI_OCN%SSI
      SNOWI => SI_OCN%SNOWI
      AUSI => ATMICE%USI
      AVSI => ATMICE%VSI
#ifdef TRACERS_WATER
      TRSI => SI_OCN%TRSI
#endif

      ATMICE%MUSI(:,:)=0
      ATMICE%HUSI(:,:)=0
      ATMICE%SUSI(:,:)=0
      ATMICE%MVSI(:,:)=0
      ATMICE%HVSI(:,:)=0
      ATMICE%SVSI(:,:)=0
#ifdef TRACERS_WATER
      ATMICE%TUSI(:,:,:)=0
      ATMICE%TVSI(:,:,:)=0
#endif

C**** Regularise ice concentration gradients to prevent advection errors
      DO J=J_0S, J_1S
      DO I=1,IM
        IF (RSI(I,J).gt.1d-4) THEN
          IF (RSISAVE(I,J).gt.RSI(I,J)) THEN ! reduce gradients
            FRSI=(RSISAVE(I,J)-RSI(I,J))/RSISAVE(I,J)
            RSIX(I,J)=RSIX(I,J)*(1.-FRSI)
            RSIY(I,J)=RSIY(I,J)*(1.-FRSI)
          END IF
          IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)

          IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
        ELSE
          RSIX(I,J) = 0.  ; RSIY(I,J) = 0.
        END IF
      END DO
      END DO
C**** set up local MHS array to contain all advected quantities
C**** Currently this is on atmospheric grid
      SNOWL=0.; HSNOW=0.; HICE=0.; SICEg=0.; TSNW=0.; TSIL=0.; MICE=0.

      DO J=J_0, J_1
        DO I=1,IM
      call get_snow_ice_layer(SNOWI(I,J),MSI(I,J),HSI(:,I,J),SSI(:,I,J),
#ifdef TRACERS_WATER
     *     TRSI(:,:,I,J),TRSNOW,TRICE,
#endif
     *     SNOWL,HSNOW,HICE,SICEg,TSNW,TSIL,MICE,.false.)

C-- MASS: ICE(LMI)
          MHS(1      :LMI    ,I,J) = MICE
C-- MASS: SNOW(2)
          MHS(LMI+1  :LMI+2  ,I,J) = SNOWL
C-- HEAT: ICE(LMI)
          MHS(3+LMI  :2+2*LMI,I,J) = HICE
C-- HEAT: SNOW(2)
          MHS(3+2*LMI:4+2*LMI,I,J) = HSNOW
C-- SALT: ICE(LMI)
          MHS(5+2*LMI:4+3*LMI,I,J) = SICEg
C-- SALT: SNOW (2, always =0)
          MHS(5+3*LMI:6+3*LMI,I,J) = 0.

#ifdef TRACERS_WATER
C**** add tracers to advected arrays
          DO ITR=1,NTM
            MHS(1+(3+ITR-1)*(LMI+2)    :LMI+(3+ITR-1)*(LMI+2),I,J)
     *        = TRICE(ITR,:)
            MHS(1+LMI+(3+ITR-1)*(LMI+2):(3+ITR)*(LMI+2)      ,I,J)
     *        = TRSNOW(ITR,:)
          ENDDO
#endif
        ENDDO
      ENDDO

C**** define inverse area array
      DO J=J_0, J_1
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0) THEN
          BYFOA(I,J)=BYDXYP(J)/FOCEAN(I,J)
          ATMICE%HSICNV(I,J) = RSI(I,J) * Sum(HSI(:,I,J))
        ELSE
          BYFOA(I,J)=0.
        END IF
      END DO
      END DO
C****
C**** North-South Advection of Sea Ice
C****
      SFASI  = 0.
      SFMSI(1:NTRICE) = 0.

C**** Update halo of RSI,FOCEAN
      CALL HALO_UPDATE(si_ocn%grid, RSI  , FROM=NORTH)
      CALL HALO_UPDATE(si_ocn%grid, RSISAVE  , FROM=NORTH)
      CALL HALO_UPDATE(si_ocn%grid,FOCEAN, FROM=NORTH+SOUTH)
#ifdef EXPEL_COASTAL_ICEXS
      CALL HALO_UPDATE(si_ocn%grid, MSI  , FROM=NORTH)
#endif

C**** calculate mass fluxes for the ice advection
      DO J=J_0,J_1S
#ifdef EXPEL_COASTAL_ICEXS
        if(j.ne.1 .and. j.ne.jm) then
          coastfacx =
     &         1d-3        ! kg/m2 ice mass -> ice thickness
     &        *1d-1        ! 10 cm/s speed for 1 m thickness difference
     &        *1d5/dxp(j)  ! over 100 km
          coastfacy =
     &         1d-3
     &        *1d-1
     &        *1d5/dyp(j)
        else
          coastfacx = 0.
          coastfacy = 0.
        endif
#endif
        I=IM
        DO IP1=1,IM
          USIDT(I,J)=0.
          IF (FOCEAN(I,J).gt.0 .and. FOCEAN(IP1,J).gt.0. .and.
     &         RSISAVE(I,J)+RSISAVE(IP1,J).gt.1d-4) then
            USIDT(I,J)=0.5*(ausi(i,j-1)+ausi(i,j))*dts
#ifdef EXPEL_COASTAL_ICEXS
            if(connect(i,j)+connect(ip1,j) .lt. 30) then
              cxi   = mod(int(connect(i,  j)),2)   ! nbr w of i  ?
              cxip1 = mod(int(connect(ip1,j)),4)/2 ! nbr e of ip1?
              du = (msi(i,j)-msi(ip1,j))*coastfacx
              du = min(10d0,max(-10d0,du))
              if(cxi.lt.cxip1) then
                du = max(0d0,du)
              elseif(cxi.gt.cxip1) then
                du = min(0d0,du)
              endif
              usidt(i,j) = usidt(i,j) + dts*du
            endif
#endif
          ENDIF
          I=IP1
        END DO
        IM1=IM
        DO I=1,IM
          VSIDT(I,J)=0.
          IF (FOCEAN(I,J+1).gt.0 .and. FOCEAN(I,J).gt.0. .and.
     &         RSISAVE(I,J)+RSISAVE(I,J+1).gt.1d-4) THEN
            VSIDT(I,J)=0.5*(avsi(im1,j)+avsi(i,j))*dts
#ifdef EXPEL_COASTAL_ICEXS
            if(connect(i,j)+connect(i,j+1) .lt. 30) then
              cyj   = mod(int(connect(i,j))/4,2) ! nbr s of j  ?
              cyjp1 = connect(i,j+1)/8           ! nbr n of j+1?
              dv = (msi(i,j)-msi(i,j+1))*coastfacy
              dv = min(10d0,max(-10d0,dv))
              if(cyj.lt.cyjp1) then
                dv = max(0d0,dv)
              elseif(cyj.gt.cyjp1) then
                dv = min(0d0,dv)
              endif
              vsidt(i,j) = vsidt(i,j) + dts*dv
            endif
#endif
          ENDIF
          IM1=I
        END DO
      END DO
      IF (hasNorthPole(si_ocn%grid)) THEN
        VSIDT(1:IM,JM)=0.
        USIDT(1:IM,JM)=AUSI(1,JM)*DTS
      END IF

C**** update RSISAVE for diagnostics
      RSISAVE(:,:)=RSI(:,:)

C**** Update halo of RSIY,RSIX,and MHS
      CALL HALO_UPDATE_COLUMN(si_ocn%grid, MHS  , FROM=NORTH)

      CALL HALO_UPDATE(si_ocn%grid, VSIDT, FROM=SOUTH)

cmpi MPI tag sequences in "grid_MIC" and "grid" occasionally collide
cmpi since "grid_MIC" is initialized to "grid" (so, same communicator).
cmpi Use "grid" everywhere
cmpi      CALL HALO_UPDATE(grid_MIC, RSIY, FROM=NORTH)
cmpi      CALL HALO_UPDATE(grid_MIC, RSIX, FROM=NORTH)
      CALL HALO_UPDATE(si_ocn%grid, RSIY, FROM=NORTH)
      CALL HALO_UPDATE(si_ocn%grid, RSIX, FROM=NORTH)

      DO I=1,IM
C****
C**** Calculate south-north sea ice fluxes at grid box edges
C****
      DO 120 J=J_0S,MIN(JM-2,J_1)
      IF(VSIDT(I,J).eq.0.)  GO TO 120
      FAW(I,J) = VSIDT(I,J)*DXV(J+1) ! be careful with atm. grid index
      IF(VSIDT(I,J).le.0.) THEN
C**** Sea ice velocity is southward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(I,J+1)-
     *       (1d0+FAW(I,J)*BYDXYP(J+1))*RSIY(I,J+1))*FOCEAN(I,J+1)
        FXSI(I,J)=FAW(I,J)*RSIX(I,J+1)*FOCEAN(I,J+1)
        FYSI(I,J)=FAW(I,J)*(FAW(I,J)*BYDXYP(J+1)*
     *       FAW(I,J)*RSIY(I,J+1)*FOCEAN(I,J+1) - 3d0*FASI(I,J))
        FMSJ(I,1:NTRICE,J) = FASI(I,J)*MHS(1:NTRICE,I,J+1)
      ELSE
C**** Sea ice velocity is northward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(I,J)+(1d0-FAW(I,J)*BYDXYP(J))*RSIY(I,J))
     *       *FOCEAN(I,J)
        FXSI(I,J)=FAW(I,J)*RSIX(I,J)*FOCEAN(I,J)
        FYSI(I,J)=FAW(I,J)*
     *       (FAW(I,J)*BYDXYP(J)*FAW(I,J)*RSIY(I,J)*FOCEAN(I,J)
     *       -3d0*FASI(I,J))
        FMSJ(I,1:NTRICE,J) = FASI(I,J)*MHS(1:NTRICE,I,J)
      END IF
        ATMICE%MVSI(I,J)=SUM(FMSJ(I,1:LMI+2,J))
        ATMICE%HVSI(I,J)=SUM(FMSJ(I,3+LMI:4+2*LMI,J))
        ATMICE%SVSI(I,J)=SUM(FMSJ(I,5+2*LMI:6+3*LMI,J))
#ifdef TRACERS_WATER
         DO ITR=1,NTM
           ATMICE%TVSI(I,J,ITR)=
     &          SUM(FMSJ(I,1+(3+ITR-1)*(LMI+2):(3+ITR)*(LMI+2),J))
         END DO
#endif
  120 CONTINUE
C****
C**** Calculate south-north sea ice fluxes near North Pole
C****
      IF (HAVE_NORTH_POLE) THEN
        IF(VSIDT(I,JM-1).eq.0.) cycle
        FAW(I,JM-1) = VSIDT(I,JM-1)*DXV(JM) ! careful with atm.grid index!
        IF(VSIDT(I,JM-1).le.0.) THEN
C**** Sea ice velocity is southward from North Pole box
          FASI(I,JM-1) = FAW(I,JM-1)*RSI(1,JM)*FOCEAN(1,JM)
          FXSI(I,JM-1) = 0.
          FYSI(I,JM-1) = -FAW(I,JM-1)*FASI(I,JM-1)
          FMSJ(I,1:NTRICE,JM-1) = FASI(I,JM-1)*MHS(1:NTRICE,1,JM)
        ELSE
C**** Sea ice velocity is northward into North Pole box
          FASI(I,JM-1) = FAW(I,JM-1)*FOCEAN(I,JM-1)*
     *       (RSI(I,JM-1)+(1d0-FAW(I,JM-1)*BYDXYP(JM-1))*RSIY(I,JM-1))
          FXSI(I,JM-1) = FAW(I,JM-1)*RSIX(I,JM-1)*FOCEAN(I,JM-1)
          FYSI(I,JM-1) = FAW(I,JM-1)*(FAW(I,JM-1)*BYDXYP(JM-1)*
     *       FAW(I,JM-1)*RSIY(I,JM-1)*FOCEAN(I,JM-1)-3d0*FASI(I,JM-1))
          FMSJ(I,1:NTRICE,JM-1) = FASI(I,JM-1)*MHS(1:NTRICE,I,JM-1)
        END IF
C**** Accumulate sea ice leaving and entering North Pole box
        SFASI = SFASI + FASI(I,JM-1)
        SFMSI(1:NTRICE) = SFMSI(1:NTRICE) + FMSJ(I,1:NTRICE,JM-1)
         ATMICE%MVSI(I,JM-1)=SUM(FMSJ(I,1:LMI+2,JM-1))
         ATMICE%HVSI(I,JM-1)=SUM(FMSJ(I,3+LMI:4+2*LMI,JM-1))
         ATMICE%SVSI(I,JM-1)=SUM(FMSJ(I,5+2*LMI:6+3*LMI,JM-1))
#ifdef TRACERS_WATER
           DO ITR=1,NTM
             ATMICE%TVSI(I,JM-1,ITR)=
     &            SUM(FMSJ(I,1+(3+ITR-1)*(LMI+2):(3+ITR)*(LMI+2),JM-1))
           END DO
#endif
      ENDIF

      END DO ! I loop
C****
C**** Update sea ice variables due to south-north fluxes
C****

C****Update halo of VSIDT, FASI, FAW, FOCEAN, FMSI, and FXSI
      CALL HALO_UPDATE_COLUMN(si_ocn%grid, FMSJ, FROM=SOUTH)
      CALL HALO_UPDATE(si_ocn%grid, FASI, FROM=SOUTH)
      CALL HALO_UPDATE(si_ocn%grid, FXSI, FROM=SOUTH)
      CALL HALO_UPDATE(si_ocn%grid, FYSI, FROM=SOUTH)
      CALL HALO_UPDATE(si_ocn%grid, FAW,  FROM=SOUTH)

      DO I=1,IM
      DO 330 J=J_0S, J_1S
      IF(VSIDT(I,J-1)) 240,210,280
C**** VSIDT(J-1)=0.
  210 IF(VSIDT(I,J))  220,330,230
C**** VSIDT(J-1)=0, VSIDT(J)<0.
  220 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) -  FASI(I,J)
      DO 225 K=1,NTRICE
  225 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) - FMSJ(I,K,J)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) - FYSI(I,J)
     *  + 3d0*(FAW(I,J)*ASI-DXYP(J)*FASI(I,J))) / (DXYP(J)-FAW(I,J))
      RSI(I,J)  = ASI*BYFOA(I,J)
      IF( RSI(I,J) > 1.d0  ) RSI(I,J) = 1.d0 ! avoid to round off
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) - FXSI(I,J)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 310
C**** VSIDT(J-1)=0, VSIDT(J)>0.
  230 RSI(I,J)  =  RSI(I,J) -  FASI(I,J)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0-FAW(I,J)*BYDXYP(J))
      RSIY(I,J) = RSIY(I,J)*(1d0-FAW(I,J)*BYDXYP(J))**2
      GO TO 310
C**** VSIDT(J-1)<0.
  240 IF(VSIDT(I,J))  260,250,270
C**** VSIDT(J-1)<0, VSIDT(J)=0.
  250 RSI(I,J)  =  RSI(I,J) +  FASI(I,J-1)*BYFOA(I,J)
      TMP = (1d0+FAW(I,J-1)*FOCEAN(I,J-1)*BYFOA(I,J))
      RSIX(I,J) = RSIX(I,J)*TMP
      RSIY(I,J) = RSIY(I,J)*TMP**2
      GO TO 310
C**** VSIDT(J-1)<0, VSIDT(J)<0  or  VSIDT(J-1)>0, VSIDT(J) not 0.
  260 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + ( FASI(I,J-1)- FASI(I,J))
      DO 265 K=1,NTRICE
  265 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) +
     *       (FMSJ(I,K,J-1)-FMSJ(I,K,J))
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J)+
     *     (FYSI(I,J-1)-FYSI(I,J)) + 3d0*((FAW(I,J-1)+
     *     FAW(I,J))*ASI-DXYP(J)*(FASI(I,J-1)+FASI(I,J))))
     *    / (DXYP(J) + (FAW(I,J-1)-FAW(I,J)))
      RSI(I,J)  = ASI*BYFOA(I,J)
      IF( RSI(I,J) > 1.d0  ) RSI(I,J) = 1.d0 ! avoid to round off
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) + (FXSI(I,J-1)-FXSI(I,J))*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 310
C**** VSIDT(J-1)<0, VSIDT(J)>0.
  270 RSI(I,J)  =  RSI(I,J) + (FASI(I,J-1)-FASI(I,J))*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+(FAW(I,J-1)*FOCEAN(I,J-1)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))
      RSIY(I,J) = RSIY(I,J)*(1d0+(FAW(I,J-1)*FOCEAN(I,J-1)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))**2
      GO TO 310
C**** VSIDT(J-1)>0.
  280 IF(VSIDT(I,J).ne.0.)  GO TO 260
C**** VSIDT(J-1)>0, VSIDT(J)=0.
      ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + FASI(I,J-1)
      DO 285 K=1,NTRICE
  285 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) + FMSJ(I,K,J-1)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) + FYSI(I,J-1)
     *    + 3d0*(FAW(I,J-1)*ASI-DXYP(J)*FASI(I,J-1))) /
     *     (DXYP(J)+FAW(I,J-1))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) + FXSI(I,J-1)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
C**** Limit RSIX and RSIY so that sea ice is positive at the edges
  310 RSI(I,J) = MAX(0d0,RSI(I,J))
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
      GO TO 330

C**** Sea ice crunches into itself and completely covers grid box
C**** Move excess ice to lower layers, allow snow to pile up
  320 RSI(I,J)   = 1d0
      RSIX(I,J)  = 0.
      RSIY(I,J)  = 0.
      DO K=1,NTRICE/(LMI+2)
C**** K=1: MASS, K=2: HEAT, K=3: SALT
        MHS(1+(LMI+2)*(K-1),I,J) = AMSI(1+(LMI+2)*(K-1)) / ASI          # 4 ice layers
        MHS(2+(LMI+2)*(K-1),I,J) = AMSI(2+(LMI+2)*(K-1)) / ASI
        DMHSI = SUM(AMSI(1+(LMI+2)*(K-1):4+(LMI+2)*(K-1)))
     *            *(BYFOA(I,J) -1d0/ ASI)
        MHS(3+(LMI+2)*(K-1),I,J) = AMSI(3+(LMI+2)*(K-1)) / ASI +
     *        XSI(3)*DMHSI
        MHS(4+(LMI+2)*(K-1),I,J) = AMSI(4+(LMI+2)*(K-1)) / ASI +
     *        XSI(4)*DMHSI
        MHS(5+(LMI+2)*(K-1),I,J) = AMSI(5+(LMI+2)*(K-1)) * BYFOA(I,J)   # 2 snow layers
        MHS(6+(LMI+2)*(K-1),I,J) = AMSI(6+(LMI+2)*(K-1)) * BYFOA(I,J)
      END DO

C**** End of loop over J
  330 CONTINUE
C**** End of loop over I
      END DO
C****
C**** Advection of Sea Ice leaving and entering North Pole box
C****
      IF (HAVE_NORTH_POLE) THEN
        ASI = RSI(1,JM)*DXYP(JM)*FOCEAN(1,JM) + SFASI/IM
        DO 345 K=1,NTRICE
  345   AMSI(K) = RSI(1,JM)*DXYP(JM)*MHS(K,1,JM)*FOCEAN(1,JM)+
     &            SFMSI(K)/IM
        IF(ASI.gt.DXYP(JM)*FOCEAN(1,JM))  GO TO 350
        RSI(1,JM)   = ASI*BYFOA(1,JM)
        IF (ASI.gt.0) MHS(1:NTRICE,1,JM) = AMSI(1:NTRICE)/ASI
        GO TO 400

C**** Sea ice crunches into itself at North Pole box
C**** Move excess ice to lower layers, allow snow to pile up
  350 RSI(1,JM)   = 1d0
      DO K=1,NTRICE/(LMI+2)
C**** K=1: MASS, K=2: HEAT, K=3: SALT
        MHS(1+(LMI+2)*(K-1),1,JM) = AMSI(1+(LMI+2)*(K-1)) / ASI         # 4 ice layers
        MHS(2+(LMI+2)*(K-1),1,JM) = AMSI(2+(LMI+2)*(K-1)) / ASI
        DMHSI = SUM(AMSI(1+(LMI+2)*(K-1):4+(LMI+2)*(K-1)))
     *            *(BYFOA(1,JM) -1d0/ ASI)
        MHS(3+(LMI+2)*(K-1),1,JM) = AMSI(3+(LMI+2)*(K-1)) / ASI +
     *        XSI(3)*DMHSI
        MHS(4+(LMI+2)*(K-1),1,JM) = AMSI(4+(LMI+2)*(K-1)) / ASI +
     *        XSI(4)*DMHSI
        MHS(5+(LMI+2)*(K-1),1,JM) = AMSI(5+(LMI+2)*(K-1)) * BYFOA(1,JM) # 2 snow layers
        MHS(6+(LMI+2)*(K-1),1,JM) = AMSI(6+(LMI+2)*(K-1)) * BYFOA(1,JM)
      END DO

      END IF   !HAVE_NORTH_POLE
C****
C**** East-West Advection of Sea Ice
C****
  400 DO 640 J=J_0S, J_1S
C****
C**** Calculate west-east sea ice fluxes at grid box edges
C****
      I=IM
      DO IP1=1,IM
      IF(USIDT(I,J).eq.0.)  GO TO 420
      FAW(I,J) = USIDT(I,J)*DYP(J)
      IF(USIDT(I,J).le.0.) THEN
C**** Sea ice velocity is westward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(IP1,J)-
     *       (1d0+FAW(I,J)*BYDXYP(J))*RSIX(IP1,J))
     *       *FOCEAN(IP1,J)
        FXSI(I,J)=FAW(I,J)*(FAW(I,J)*BYDXYP(J)*
     *       FAW(I,J)*RSIX(IP1,J)*FOCEAN(IP1,J)-3d0*FASI(I,J))
        FYSI(I,J)=FAW(I,J)*RSIY(IP1,J)*FOCEAN(IP1,J)
        FMSI(1:NTRICE,I) = FASI(I,J)*MHS(1:NTRICE,IP1,J)
      ELSE
C**** Sea ice velocity is eastward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(I,J)+(1d0-FAW(I,J)*BYDXYP(J))*RSIX(I,J))
     *       *FOCEAN(I,J)
        FXSI(I,J)=FAW(I,J)*(FAW(I,J)*BYDXYP(J)*FAW(I,J)*
     *       RSIX(I,J)*FOCEAN(I,J)-3d0*FASI(I,J))
        FYSI(I,J)=FAW(I,J)*RSIY(I,J)*FOCEAN(I,J)
        FMSI(1:NTRICE,I) = FASI(I,J)*MHS(1:NTRICE,I,J)
      END IF
         ATMICE%MUSI(I,J)=SUM(FMSI(1:LMI+2,I))
         ATMICE%HUSI(I,J)=SUM(FMSI(3+LMI:4+2*LMI,I))
         ATMICE%SUSI(I,J)=SUM(FMSI(5+2*LMI:6+3*LMI,I))
#ifdef TRACERS_WATER
         DO ITR=1,NTM
           ATMICE%TUSI(I,J,ITR)=
     &          SUM(FMSI(1+(3+ITR-1)*(LMI+2):(3+ITR)*(LMI+2),I))
         END DO
#endif
  420 I=IP1
      END DO
C****
C**** Update sea ice variables due to west-east fluxes
C****
      IM1=IM
      DO 630 I=1,IM
      IF(USIDT(IM1,J)) 540,510,580
C**** USIDT(IM1)=0.
  510 IF(USIDT(I,J))  520,630,530
C**** USIDT(IM1)=0, USIDT(I)<0.
  520 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) -  FASI(I,J)
      DO 525 K=1,NTRICE
  525 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) - FMSI(K,I)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) - FXSI(I,J)
     *    + 3d0*(FAW(I,J)*ASI-DXYP(J)*FASI(I,J))) / (DXYP(J)-FAW(I,J))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) - FYSI(I,J)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 610
C**** USIDT(IM1)=0, USIDT(I)>0.
  530 RSI(I,J)  =  RSI(I,J) -  FASI(I,J)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0-FAW(I,J)*BYDXYP(J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0-FAW(I,J)*BYDXYP(J))
      GO TO 610
C**** USIDT(IM1)<0.
  540 IF(USIDT(I,J))  560,550,570
C**** USIDT(IM1)<0, USIDT(I)=0.
  550 RSI(I,J)  =  RSI(I,J) +  FASI(IM1,J)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+FAW(IM1,J)*FOCEAN(IM1,J)*BYFOA(I,J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0+FAW(IM1,J)*FOCEAN(IM1,J)*BYFOA(I,J))
      GO TO 610
C**** USIDT(IM1)<0, USIDT(I)<0  or  USIDT(IM1)>0, USIDT(I) not 0.
  560 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + (FASI(IM1,J)- FASI(I,J))
      DO 565 K=1,NTRICE
  565 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) +
     *       (FMSI(K,IM1)-FMSI(K,I))
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J)+
     *     (FXSI(IM1,J)-FXSI(I,J)) + 3d0*((FAW(IM1,J)+FAW(I,J))*ASI-
     *     DXYP(J)*(FASI(IM1,J)+FASI(I,J))))
     *    / (DXYP(J) + (FAW(IM1,J)-FAW(I,J)))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) + (FYSI(IM1,J)-FYSI(I,J))*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 610
C**** USIDT(IM1)<0, USIDT(I)>0.
  570 RSI(I,J)  =  RSI(I,J) + (FASI(IM1,J)-FASI(I,J))*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+(FAW(IM1,J)*FOCEAN(IM1,J)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0+(FAW(IM1,J)*FOCEAN(IM1,J)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))
      GO TO 610
C**** USIDT(IM1)>0.
  580 IF(USIDT(I,J).ne.0.)  GO TO 560
C**** USIDT(IM1)>0, USIDT(I)=0.
      ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + FASI(IM1,J)
      DO 585 K=1,NTRICE
  585 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) + FMSI(K,IM1)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) + FXSI(IM1,J)
     *    + 3d0*(FAW(IM1,J)*ASI-DXYP(J)*FASI(IM1,J))) /
     *     (DXYP(J)+FAW(IM1,J))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) + FYSI(IM1,J)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
C**** Limit RSIX and RSIY so that sea ice is positive at the edges
  610 RSI(I,J) = MAX(0d0,RSI(I,J))
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)>1d0) RSI(I,J)=1d0
      GO TO 630

C**** Sea ice crunches into itself and completely covers grid box
C**** Move excess ice to lower layers, allow snow to pile up
  620 RSI(I,J)   = 1d0
      RSIX(I,J)  = 0.
      RSIY(I,J)  = 0.
      DO K=1,NTRICE/(LMI+2)
C**** K=1: MASS, K=2: HEAT, K=3: SALT
        MHS(1+(LMI+2)*(K-1),I,J) = AMSI(1+(LMI+2)*(K-1)) / ASI          # 4 ice layers
        MHS(2+(LMI+2)*(K-1),I,J) = AMSI(2+(LMI+2)*(K-1)) / ASI
        DMHSI = SUM(AMSI(1+(LMI+2)*(K-1):4+(LMI+2)*(K-1)))
     *            *(BYFOA(I,J) -1d0/ ASI)
        MHS(3+(LMI+2)*(K-1),I,J) = AMSI(3+(LMI+2)*(K-1)) / ASI +
     *        XSI(3)*DMHSI
        MHS(4+(LMI+2)*(K-1),I,J) = AMSI(4+(LMI+2)*(K-1)) / ASI +
     *        XSI(4)*DMHSI
        MHS(5+(LMI+2)*(K-1),I,J) = AMSI(5+(LMI+2)*(K-1)) * BYFOA(I,J)   # 2 snow layers
        MHS(6+(LMI+2)*(K-1),I,J) = AMSI(6+(LMI+2)*(K-1)) * BYFOA(I,J)
      END DO

C**** End of loop over I
  630 IM1=I
C**** End of loop over J
  640 CONTINUE

      IF (KOCEAN.ge.1) THEN ! full ocean calculation, adjust sea ice
C**** set global variables from local array
C**** Currently on atmospheric grid, so no interpolation necessary
        DO J=J_0, J_1
          DO I=1,si_ocn%IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
C**** Fresh water sea ice mass convergence (needed for qflux model)
            atmice%MSICNV(I,J) = RSI(I,J)*(SUM(MHS(1:LMI+2,I,J))
     *           -SUM(MHS(5+2*LMI:4+3*LMI,I,J)))
     *                  - RSISAVE(I,J)*(ACE1I+SNOWI(I,J)
     *           +MSI(I,J)-SUM(SSI(1:LMI,I,J)))

C-- MASS: ICE(LMI)
          MICE  = MHS(1      :LMI    ,I,J)
C-- MASS: SNOW(2)
          SNOWL = MHS(LMI+1  :LMI+2  ,I,J)
C-- HEAT: ICE(LMI)
          HICE  = MHS(3+LMI  :2+2*LMI,I,J)
C-- HEAT: SNOW(2)
          HSNOW = MHS(3+2*LMI:4+2*LMI,I,J)
C-- SALT: ICE(LMI)
          SICEg = MHS(5+2*LMI:4+3*LMI,I,J)
C-- SALT: SNOW (2, always =0)
c         MHS(5+3*LMI:6+3*LMI,I,J) = 0.

#ifdef TRACERS_WATER
          DO ITR=1,NTM
            TRICE(ITR,:) =
     *        MHS(1+(3+ITR-1)*(LMI+2)    :LMI+(3+ITR-1)*(LMI+2),I,J)
            TRSNOW(ITR,:) =
     *        MHS(1+LMI+(3+ITR-1)*(LMI+2):(3+ITR)*(LMI+2)      ,I,J)
          ENDDO
#endif

C-- Mass fluxes required to keep first layer ice = ACE1I
      FMSI2=MICE(1)+MICE(2)-ACE1I
C**** relayer lower ice layers
        call relayer(FMSI2,MICE,HICE,SICEg
#ifdef TRACERS_WATER
     *       ,TRICE
#endif
     *       )

C**** relayer upper two layers
        call relayer_12(HSNOW,HICE,SICEg,MICE,SNOWL
#ifdef TRACERS_WATER
     *       ,TRSNOW,TRICE
#endif
     *       )

C**** reconstitute upper snow and ice layers
        call set_snow_ice_layer(HSNOW,HICE,SICEg,MICE,SNOWL,
#ifdef TRACERS_WATER
     *       TRSNOW,TRICE,TRSI(:,:,I,J),
#endif
     *       SNOWI(I,J),MSI1,MSI(I,J),HSI(:,I,J),SSI(:,I,J))

            ATMICE%HSICNV(I,J) = FOCEAN(I,J) *
     *         (RSI(I,J)*Sum(HSI(:,I,J)) - ATMICE%HSICNV(I,J))

            atmice%FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)-
     *           SUM(SSI(1:LMI,I,J)))
            END IF
          END DO
        END DO
cC**** Set atmospheric arrays
c        DO J=J_0, J_1
c          DO I=1,si_ocn%IMAXJ(J)
c            IF (FOCEAN(I,J).gt.0) THEN
cC**** adjust rad fluxes for change in ice fraction
c              if (rsi(i,j).gt.rsisave(i,j)) ! ice from ocean
c     *       call RESET_SURF_FLUXES(I,J,1,2,RSISAVE(I,J),RSI(I,J))
c              if (rsi(i,j).lt.rsisave(i,j)) ! ocean from ice
c     *       call RESET_SURF_FLUXES(I,J,2,1,1.-RSISAVE(I,J),1.-RSI(I,J))
cC****
c            END IF
c          END DO
c        END DO
      ELSE          ! fixed SST case, save implied heat convergence
        DO J=J_0, J_1
          DO I=1,si_ocn%IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
              atmice%HSICNV(I,J)=RSI(I,J)*SUM(MHS(3+LMI:4+2*LMI,I,J)) -
     -                              ATMICE%HSICNV(I,J)
C**** reset sea ice concentration
              RSI(I,J)=RSISAVE(I,J)
            END IF
          END DO
        END DO
      END IF
C****

      DEALLOCATE(FMSI,FMSJ,MHS,SFMSI,AMSI)
      DEALLOCATE(FASI, FXSI, FYSI, FAW, BYFOA, USIDT, VSIDT)

      RETURN
      END SUBROUTINE ADVSI
#endif

#ifdef CUBED_SPHERE
      subroutine INT_AtmA2IceA_XY(aA,iA)
!@sum interpolate from Atm A-grid to Ice A-grid
!@+   CS version uses interpolation in CS XY-space
!@auth Denis Gueyffier
!@auth M. Kelley (cs2llint_ij)
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : ICE_PACK=>PACK_DATA,am_i_root,
     &     HALO_UPDATE
      USE ICEDYN, only : grid_ICDYN,IMICDYN,JMICDYN
      USE ICEDYN_COM, only : CS2ICEint_a
      USE cs2ll_utils, only : cs2llint_ij
      IMPLICIT NONE
      real*8 ::
     &     aA(agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO),
     &     iA(1:IMICDYN,
     &     grid_ICDYN%J_STRT_HALO:grid_ICDYN%J_STOP_HALO)
c      character*80 :: title
c      real*8, allocatable :: iA_glob(:,:)
c      real*4, allocatable :: iA4_glob(:,:)
c      allocate(iA_glob(IMICDYN,JMICDYN),iA4_glob(IMICDYN,JMICDYN))
      call cs2llint_ij(agrid,CS2ICEint_a,aA,iA)
c      call HALO_UPDATE(grid_ICDYN,iA)
c      call ICE_PACK(grid_ICDYN,iA,iA_glob)
c      if (am_i_root()) then
c      open(900,FILE="aA2iA",FORM='unformatted',
c     &        STATUS='unknown')
c      title="test"
c      iA4_glob=iA_glob
c      write(900) title,iA4_glob
c      close(900)
c      endif
c      deallocate(iA_glob,iA4_glob)

      end subroutine INT_AtmA2IceA_XY
#else
c***  for the moment me assume that the atm and icedyn grids are
c***  both latlon with equal resolution
c      iA=aA
#endif

      SUBROUTINE init_icedyn(iniOCEAN,atmice)
!@sum  init_icedyn initializes ice dynamics variables
!@auth Gavin Schmidt
      USE MODEL_COM, only : dtsrc
      USE MDIAG_COM, only : ia_src=>ia_cpl
      USE MDIAG_COM, only : make_timeaxis
      USE DOMAIN_DECOMP_1D, only : getDomainBounds,ICE_HALO=>HALO_UPDATE
     &     ,hasSouthPole,hasNorthPole
      USE ICEDYN_COM, only : igice
     &     ,kicij,ia_icij,denom_icij,igrid_icij,jgrid_icij,lname_icij
     &     ,sname_icij,units_icij,scale_icij,ij_usi,ij_vsi,ij_dmui
     &     ,ij_dmvi,ij_rsi,ij_pice,ij_sispeed
     &     ,cdl_icij
      USE ICEDYN, only : ifocean=>focean,
     &     osurf_tilt,bydts,usi,vsi,uice,vice,lon_dg,lat_dg
      USE ICEDYN, only : NX1,grid_ICDYN,grid_NXY,IMICDYN,JMICDYN,
     &     GEOMICDYN,ICDYN_MASKS,DLATM
#ifdef CUBED_SPHERE
      use dd2d_utils, only : ATM_HALO=>HALO_UPDATE
      USE ICEDYN, only : lon,lat,lonb,latb,uvm
      USE ICEDYN_COM, only : CS2ICEint_a,CS2ICEint_b,i2a_uc,i2a_vc
     &     ,ICE2CSint ,UVLLATUC,UVLLATVC,UVMATUC,UVMATVC,CONNECT
     &     ,FOA,BYFOA
      USE cs2ll_utils, only : init_cs2llint_type,init_ll2csint_type,
     &     ll2csint_ij
      USE GEOM, only : AXYP,BYAXYP,lon2d,lat2d,lonuc,latuc,lonvc,latvc
      use constant, only : pi
c      USE DOMAIN_DECOMP_1D, only : READT_PARALLEL
c      USE FILEMANAGER, only : openunit,closeunit,nameunit
#else
      use domain_decomp_1d, only : init_band_pack_type,band_pack
      use icedyn_com, only : pack_a2i,pack_i2a
#ifdef EXPEL_COASTAL_ICEXS
      use domain_decomp_1d, only : ATM_HALO=>HALO_UPDATE
      use icedyn_com, only : connect
#endif
#endif
      use Dictionary_mod
      use cdl_mod
      USE EXCHANGE_TYPES, only : atmice_xchng_vars
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: iniOCEAN
      type(atmice_xchng_vars) :: atmice
c
      INTEGER i,j,k,kk,J_0,J_1,J_0H,J_1H,J_1S,im1,ip1,im
      character(len=10) :: xstr,ystr
      logical :: set_miss
#ifdef CUBED_SPHERE
      integer :: imin,imax,jmin,jmax,iu_mask
      real*8 :: lonb_tmp(imicdyn)
      real*8, dimension(:,:), allocatable :: uvm_tmp
#endif

C**** First, set up the ice dynamics lat-lon grid.
C**** The resolutions IMICDYN, JMICDYN are defined in ICEDYN.f.

C**** Calculate spherical geometry
      call GEOMICDYN()
      igice%dlatm = dlatm

#ifdef CUBED_SPHERE
c Get the ice dynamics land mask from the ocean topo file
c      call openunit("TOPO_OC",iu_mask,.true.,.true.)
c      CALL READT_PARALLEL(grid_icdyn,iu_mask,NAMEUNIT(iu_mask),
c     &     iFOCEAN,1)
c      call closeunit(iu_mask)
c**** set up CS2ICEint, a data structure for CS to latlon interpolation
      call init_cs2llint_type(atmice%grid,grid_ICDYN,lon,lat,1,JMICDYN,
     &     CS2ICEint_a)
      lonb_tmp = lonb-pi
      call init_cs2llint_type(atmice%grid,grid_ICDYN,lonb_tmp,latb,
     &     1,JMICDYN-1,
     &     CS2ICEint_b,setup_rot_pol=.true.)
C**** Derive the ice dynamics land mask from the atmosphere mask
      call INT_AtmA2IceA_XY(atmice%focean,iFocean)
#else
C**** Initialize derived types used for passing arrays between the
C**** atmosphere and ice dynamics
      call init_band_pack_type(atmice%grid, grid_ICDYN,
     &     grid_ICDYN%J_STRT_HALO,grid_ICDYN%J_STOP_HALO,
     &     pack_a2i)
      call init_band_pack_type(grid_ICDYN, atmice%grid,
     &     atmice%grid%J_STRT_HALO,atmice%grid%J_STOP_HALO,
     &     pack_i2a)
      igice%pack_a2i => pack_a2i
      igice%pack_i2a => pack_i2a
C**** The ice dynamics land mask is that of the atmosphere
      call band_pack(pack_a2i, atmice%focean, ifocean) ! ifocean = afocean
#endif

      icedyn_processors_only: if(grid_ICDYN%have_domain) then
        call ICE_HALO(grid_ICDYN, iFOCEAN)
      endif icedyn_processors_only
        call ICDYN_MASKS()

      bydts = 1./dtsrc

      call getDomainBounds( grid_NXY,J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &                   J_STRT     =J_0,  J_STOP =J_1,
     &                   J_STOP_SKP =J_1S                   )

C**** Initialise ice dynamics if ocean model needs initialising
      if (iniOCEAN) THEN
         USI=0.
         VSI=0.
      endif

#ifdef CUBED_SPHERE
C**** precompute some arrays for ice advection on the atm grid
      do j=atmice%grid%j_strt,atmice%grid%j_stop
        do i=atmice%grid%i_strt,atmice%grid%i_stop
          FOA(I,J)=AXYP(I,J)*atmice%focean(I,J)
          IF(atmice%focean(I,J).gt.0) THEN
            BYFOA(I,J)=BYAXYP(I,J)/atmice%focean(I,J)
          ELSE
            BYFOA(I,J)=0.
          END IF
        enddo
      enddo
      call ATM_HALO(atmice%grid, FOA)
      call ATM_HALO(atmice%grid, BYFOA)
C**** precompute some interpolation info for ice advection
c ice b-grid -> atm "west" edges
      imin=lbound(lonuc,1); imax=ubound(lonuc,1)
      jmin=lbound(lonuc,2); jmax=ubound(lonuc,2)
      call init_ll2csint_type(grid_icdyn,atmice%grid,
     &     lonb,latb, 1,JMICDYN-1,
     &     imin,imax,jmin,jmax,lonuc,latuc,
     &     i2a_uc)
      allocate(uvllatuc(2,imin:imax,jmin:jmax))
      allocate(uvmatuc(imin:imax,jmin:jmax))
c ice b-grid -> atm "south" edges
      imin=lbound(lonvc,1); imax=ubound(lonvc,1)
      jmin=lbound(lonvc,2); jmax=ubound(lonvc,2)
      call init_ll2csint_type(grid_icdyn,atmice%grid,
     &     lonb,latb, 1,JMICDYN-1,
     &     imin,imax,jmin,jmax,lonvc,latvc,
     &     i2a_vc)
      allocate(uvllatvc(2,imin:imax,jmin:jmax))
      allocate(uvmatvc(imin:imax,jmin:jmax))
c ice b-grid -> atm a-grid
      call init_ll2csint_type(grid_icdyn,atmice%grid,
     &     lonb,latb, 1,JMICDYN-1,
     &     atmice%grid%isd,atmice%grid%ied,
     &          atmice%grid%jsd,atmice%grid%jed,
     &     lon2d,lat2d,
     &     ICE2CSint,skip_halos=.true.)

c Interpolate the ice dynamics velocity mask to atm cell edges, and
c encode the "ocean-connectedness" of gridpoint i,j using:
c ["west:" 1] + ["east:" 2] + ["south:" 4] + ["north:" 8].
      allocate(uvm_tmp(imicdyn,j_0h:j_1h))
      uvm_tmp(1:imicdyn,j_0:j_1) = uvm(2:imicdyn+1,j_0:j_1)
      call ll2csint_ij(grid_icdyn,i2a_uc,uvm_tmp,uvmatuc)
      call ll2csint_ij(grid_icdyn,i2a_vc,uvm_tmp,uvmatvc)
      deallocate(uvm_tmp)
      do j=atmice%grid%j_strt,atmice%grid%j_stop
        do i=atmice%grid%i_strt,atmice%grid%i_stop
          connect(i,j) = 0
c Gridcells whose ocean fraction is <10% are not ocean-connected.
c          if(atmice%focean(i,j).gt.0.1d0) then
c            if(atmice%focean(i-1,j).gt.0.1d0 .and. uvmatuc(i  ,j).gt.0.)
c     &           connect(i,j) = connect(i,j) + 1
c            if(atmice%focean(i+1,j).gt.0.1d0 .and. uvmatuc(i+1,j).gt.0.)
c     &           connect(i,j) = connect(i,j) + 2
c            if(atmice%focean(i,j-1).gt.0.1d0 .and. uvmatvc(i,j  ).gt.0.)
c     &           connect(i,j) = connect(i,j) + 4
c            if(atmice%focean(i,j+1).gt.0.1d0 .and. uvmatvc(i,j+1).gt.0.)
c     &           connect(i,j) = connect(i,j) + 8
c          endif
          if(atmice%focean(i,j).gt.0.0d0) then
            if(atmice%focean(i-1,j).gt.0.0d0)
     &           connect(i,j) = connect(i,j) + 1
            if(atmice%focean(i+1,j).gt.0.0d0)
     &           connect(i,j) = connect(i,j) + 2
            if(atmice%focean(i,j-1).gt.0.0d0)
     &           connect(i,j) = connect(i,j) + 4
            if(atmice%focean(i,j+1).gt.0.0d0)
     &           connect(i,j) = connect(i,j) + 8
          endif
        enddo
      enddo
      call atm_halo(atmice%grid, connect)
#else
#ifdef EXPEL_COASTAL_ICEXS
c encode the "ocean-connectedness" of gridpoint i,j using:
c ["west:" 1] + ["east:" 2] + ["south:" 4] + ["north:" 8].
      im = atmice%grid%i_stop
      do j=atmice%grid%j_strt_skp,atmice%grid%j_stop_skp
        im1 = im-1
        i   = im
        do ip1=1,im
          connect(i,j) = 0
          if(atmice%focean(i,j).gt.0.0d0) then
            if(atmice%focean(im1,j).gt.0.0d0)
     &           connect(i,j) = connect(i,j) + 1
            if(atmice%focean(ip1,j).gt.0.0d0)
     &           connect(i,j) = connect(i,j) + 2
            if(atmice%focean(i,j-1).gt.0.0d0)
     &           connect(i,j) = connect(i,j) + 4
            if(atmice%focean(i,j+1).gt.0.0d0)
     &           connect(i,j) = connect(i,j) + 8
          endif
          im1 = i
          i = ip1
        enddo
      enddo
      if(hasSouthPole(atmice%grid)) then
        connect(:, 1) = 0       ! assume landlocked
      endif
      if(hasNorthPole(atmice%grid)) then
        connect(:,atmice%grid%j_stop) = 1+2+4+8 ! assume ocean
      endif
      call atm_halo(atmice%grid, connect)
#endif
#endif

C**** set uisurf,visurf for atmospheric drag calculations
      call get_uisurf(usi,vsi,atmice)

C**** set properties for ICIJ diagnostics
      do k=1,kicij
        ia_icij(k)=ia_src
        denom_icij(k)=0
        igrid_icij(k)=1
        jgrid_icij(k)=1
      enddo
      k=0

      k=k+1
      IJ_RSI=k
      lname_icij(k)="Ocean ice fraction (ice dynamic grid)"
      sname_icij(k)="icij_rsi"
      units_icij(k)="%"
      scale_icij(k)=100.

      k=k+1
      IJ_USI=k
c      denom_icij(k)=IJ_RSIU ! need to add IJ_RSIU
      lname_icij(k)="Sea ice EW velocity x POICEU"
      sname_icij(k)="icij_usi"
      units_icij(k)="m/s"
      scale_icij(k)=1.
      igrid_icij(k)=2

      k=k+1
      IJ_VSI=k
c      denom_icij(k)=IJ_RSIV ! need to add IJ_RSIV
      lname_icij(k)="Sea ice NS velocity x POICEV"
      sname_icij(k)="icij_vsi"
      units_icij(k)="m/s"
      scale_icij(k)=1.
      jgrid_icij(k)=2

      k=k+1
      IJ_SISPEED=k
      lname_icij(k)="Sea ice speed x POICE"
      sname_icij(k)="icij_sispeed"
      units_icij(k)="m/s"
      scale_icij(k)=1.

      k=k+1
      IJ_DMUI=k
      lname_icij(k)="Ice-ocean EW stress"
      sname_icij(k)="icij_dmui"
      units_icij(k)="kg/m s^2"
      scale_icij(k)=1./dtsrc

      k=k+1
      IJ_DMVI=k
      lname_icij(k)="Ice-ocean NS stress"
      sname_icij(k)="icij_dmvi"
      units_icij(k)="kg/m s^2"
      scale_icij(k)=1./dtsrc

      k=k+1
      IJ_PICE=k
      denom_icij(k)=IJ_RSI
      lname_icij(k)="Sea ice internal pressure"
      sname_icij(k)="icij_psi"
      units_icij(k)="10^3 kg/m s^2"
      scale_icij(k)=1d-3

      if (k.gt.KICIJ) then
        write(6,*) "Too many ICIJ diags: increase KICIJ to at least"
     *       ,k
        call stop_model("ICIJ diagnostic error",255)
      end if

c
c Declare the dimensions and metadata of output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c
      call init_cdl_type('cdl_icij',cdl_icij)
      call add_coord(cdl_icij,'lon',imicdyn,units='degrees_east',
     &     coordvalues=lon_dg(:,1))
      call add_coord(cdl_icij,'lat',jmicdyn,units='degrees_north',
     &     coordvalues=lat_dg(:,1))
      call add_coord(cdl_icij,'lon2',imicdyn,units='degrees_east',
     &     coordvalues=lon_dg(:,2))
      call add_coord(cdl_icij,'lat2',jmicdyn,units='degrees_north',
     &     coordvalues=lat_dg(:,2))
      do k=1,kicij
        if(trim(sname_icij(k)).eq.'unused') cycle
        xstr='lon) ;'
        if(igrid_icij(k).eq.2) xstr='lon2) ;'
        ystr='(lat,'
        if(jgrid_icij(k).eq.2) ystr='(lat2,'
        set_miss = denom_icij(k).ne.0
        call add_var(cdl_icij,
     &       'float '//trim(sname_icij(k))//trim(ystr)//trim(xstr),
     &       units=trim(units_icij(k)),
     &       long_name=trim(lname_icij(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
      enddo

      RETURN
      END SUBROUTINE init_icedyn

#ifndef CUBED_SPHERE
      SUBROUTINE diag_ICEDYN
!!@sum  diag_ICEDYN prints out diagnostics for ice dynamics
!!@&    ESMF: It should only be called from a serial region.
!!@$          It is NOT parallelized
!!@auth Gavin Schmidt
!      USE CONSTANT, only : undef,teeny
!      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
!     *     ,amon0,jdate,amon,jyear
!      USE DIAG_COM, only : qdiag
!      USE MDIAG_COM, only : acc_period,
!     &     lname_strlen,sname_strlen,units_strlen
!      USE ICEDYN_COM
!      USE DIAG_SERIAL, only : focean=>FOCEAN_glob
!      USE FILEMANAGER, only : openunit
!      IMPLICIT NONE
!      REAL*8, DIMENSION(IMIC,JMIC) :: Q,ADENOM
!      INTEGER I,J,L,N,KXLB,IP1,k1,k
!      CHARACTER XLB*30
!      CHARACTER TITLE*80
!      CHARACTER(len=lname_strlen) :: lname
!      CHARACTER(len=sname_strlen) :: sname
!      CHARACTER(len=units_strlen) :: units
!      REAL*8 QJ(JM),QSUM
!      REAL*8 byiacc
!
!      IF (.not. QDIAG) RETURN
!C**** determine label to be added to all titles
!      KXLB = INDEX(XLABEL(1:11),'(')-1
!      IF(KXLB.le.0) KXLB = 10
!      XLB = ' '
!      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
!      XLB = TRIM(XLB)//" "//XLABEL(1:KXLB)
!
!C**** Open output files
!      call open_ij(trim(acc_period)//'.icij'//
!     *     XLABEL(1:LRUNID),imic,jmic)
!C****
!C**** Simple scaled ICIJ diagnostics
!C****
!      DO K=1,KICIJ
!        byiacc=1./(IDACC(IA_ICIJ(K))+teeny)
!        lname=lname_icij(k)
!        if(denom_icij(k).gt.0) then
!          adenom=icijg(1:imic,1:jmic,denom_icij(k)) * byiacc
!        else
!          adenom=1.
!        endif
!        k1 = index(lname,' x ')
!        if (k1 .gt. 0) then
!          if (index(lname,' x POICEU') .gt. 0) then
!            do j=1,jmic
!            i=imic
!            do ip1=1,imic
!              adenom(i,j)=0.5*(icijg(i,j,ij_rsi)+icijg(ip1,j,ij_rsi))
!     *             * byiacc
!              i=ip1
!            end do
!            end do
!          else if (index(lname,' x POICEV') .gt. 0) then
!            do j=1,jmic-1
!            do i=1,imic
!              adenom(i,j)=0.5*(icijg(i,j,ij_rsi)+icijg(i,j+1,ij_rsi))
!     *             * byiacc
!            end do
!            end do
!          end if
!          lname(k1:50) = ' '
!        end if
!
!        Q=UNDEF ; QJ=UNDEF ; QSUM=UNDEF
!        DO J=1,JMIC
!          DO I=1,IMIC
!            IF (ADENOM(I,J).gt.0 .and. FOCEAN(I,J).gt.0.5)
!     *           Q(I,J)=SCALE_ICIJ(K)*ICIJg(I,J,K)*byiacc/adenom(i,j)
!          END DO
!        END DO
!        Q(2:IMIC,JMIC)=Q(1,JMIC)
!        Q(2:IMIC,1)=Q(1,1)
!        TITLE=trim(LNAME)//" ("//trim(UNITS_ICIJ(K))//") "
!        TITLE(51:80)=XLB
!        CALL POUT_IJ(TITLE,SNAME_ICIJ(K),LNAME_ICIJ(K),UNITS_ICIJ(K),Q
!     *       ,QJ,QSUM,IGRID_ICIJ(K),JGRID_ICIJ(K))
!
!      END DO
!
!      call close_ij
!C****
      RETURN
      END SUBROUTINE diag_ICEDYN
#endif /* not CUBED_SPHERE */

      SUBROUTINE GET_UISURF(UICE,VICE,atmice)
!@sum calculate atmos. A grid winds from B grid
!@auth Gavin Schmidt, Denis Gueyffier
!@auth M. Kelley
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE ICEDYN, only : grid=>grid_icdyn,IMICDYN,JMICDYN
#ifdef CUBED_SPHERE
      use cs2ll_utils, only : ll2csint_lij
      use icedyn_com, only : ICE2CSint
#else
      USE DOMAIN_DECOMP_1D, only : band_pack,
     &     hasSouthPole, hasNorthPole
      use icedyn_com, only : pack_i2a
      use icedyn, only : cosiu,siniu
#endif
      USE EXCHANGE_TYPES, only : atmice_xchng_vars
      IMPLICIT NONE
      REAL*8, DIMENSION(IMICDYN,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO) :: UICE,VICE
      type(atmice_xchng_vars) :: atmice
c
      INTEGER :: J_0S,J_1S,I,J,IM1,IM,JM
      real*8, allocatable :: uvice(:,:,:),uvice_cs(:,:,:)

#ifdef CUBED_SPHERE
      allocate(uvice(2,IMICDYN,grid%J_STRT_HALO:grid%J_STOP_HALO))
      do j=max(1,grid%J_STRT),min(grid%J_STOP,JMICDYN-1)
        do i=1,IMICDYN
          uvice(1,i,j) = uice(i,j)
          uvice(2,i,j) = vice(i,j)
c this is not necessary since usi,vsi already have this
          IF (abs(uvice(1,i,j)).lt.1d-10) uvice(1,i,j)=0.
          IF (abs(uvice(2,i,j)).lt.1d-10) uvice(2,i,j)=0.
        enddo
      enddo
      allocate(uvice_cs(2,
     &     atmice%grid%I_STRT_HALO:atmice%grid%I_STOP_HALO,
     &     atmice%grid%J_STRT_HALO:atmice%grid%J_STOP_HALO))
      call ll2csint_lij(grid,ICE2CSint,uvice,uvice_cs,
     &     is_ll_vector=.true.)
      do j=atmice%grid%J_STRT,atmice%grid%J_STOP
        do i=atmice%grid%I_STRT,atmice%grid%I_STOP
          atmice%uisurf(i,j) = uvice_cs(1,i,j)
          atmice%visurf(i,j) = uvice_cs(2,i,j)
        enddo
      enddo
      deallocate(uvice_cs)
      deallocate(uvice)
#else
c**** We assume that ice grid and latlon atm grid have same resolution
      im = atmice%grid%im_world
      jm = atmice%grid%jm_world
      call band_pack(pack_i2a, uice, atmice%usi) ! fills halos
      call band_pack(pack_i2a, vice, atmice%vsi) ! fills halos
      call getDomainBounds(atmice%grid,J_STRT_SKP=J_0S,J_STOP_SKP=J_1S)
      do j=J_0S,J_1S
         im1 = IM
         do i=1,IM
c*** B->A : 4 points averaging
            atmice%uisurf(i,j) =
     &          0.25*(atmice%usi(i  ,j)+atmice%usi(i  ,j-1)+
     &                atmice%usi(im1,j)+atmice%usi(im1,j-1))
            atmice%visurf(i,j) =
     &           0.25*(atmice%vsi(i  ,j)+atmice%vsi(i  ,j-1)+
     &                 atmice%vsi(im1,j)+atmice%vsi(im1,j-1))
            im1 = i
         enddo
      enddo
c*** Poles
      if (hasSouthPole(atmice%grid)) then
         atmice%uisurf(1,1) = 0. ; atmice%visurf(1,1) = 0.
         do i=1,IM
           atmice%uisurf(1,1) = atmice%uisurf(1,1)
     &          +(atmice%vsi(i,1)*siniu(i))*2./IM
           atmice%visurf(1,1) = atmice%visurf(1,1)
     &          +(atmice%vsi(i,1)*cosiu(i))*2./IM
         end do
         atmice%uisurf(:,1)=atmice%uisurf(1,1)
         atmice%visurf(:,1)=atmice%visurf(1,1)
      endif

      if (hasNorthPole(atmice%grid)) then
         atmice%uisurf(1,JM) = 0. ; atmice%visurf(1,JM) = 0.
         do i=1,IM
           atmice%uisurf(1,JM) = atmice%uisurf(1,JM)
     &          -(atmice%vsi(i,JM-1)*siniu(i))*2./IM
           atmice%visurf(1,JM) = atmice%visurf(1,JM)
     &          +(atmice%vsi(i,JM-1)*cosiu(i))*2./IM
         end do
         atmice%uisurf(:,JM)=atmice%uisurf(1,JM)
         atmice%visurf(:,JM)=atmice%visurf(1,JM)
      endif
#endif

      RETURN
      END SUBROUTINE GET_UISURF
