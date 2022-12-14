#include "rundeck_opts.h"
      MODULE ATM_COM
!@sum  ATM_COM Main atmospheric variables
!@vers 2013/04/02
!@auth Original Development Team
      USE RESOLUTION, only : LM
#ifdef USE_ESMF
      USE ESMF, only: ESMF_Clock
#endif
      IMPLICIT NONE
      SAVE

#ifdef USE_ESMF
! different components can define their own clocks later
      Type (ESMF_CLOCK) :: atmclock
#endif

!@var LM_REQ Extra number of radiative equilibrium layers
      INTEGER, PARAMETER :: LM_REQ=3
!@var REQ_FAC/REQ_FAC_M factors for REQ layer pressures
      REAL*8, PARAMETER, DIMENSION(LM_REQ-1) ::
     *     REQ_FAC=(/ .5d0, .2d0 /)               ! edge
      REAL*8, PARAMETER, DIMENSION(LM_REQ) ::
     *     REQ_FAC_M=(/ .75d0, .35d0, .1d0 /),    ! mid-points
     *     REQ_FAC_D=(/ .5d0,  .3d0,  .2d0 /)     ! delta

!@var PL00, PMIDL00, PDSIGL00, AML00 press (mb), mid-pressure (mb),
!@+        pressure thickness (mb), mass (kg/m2) for mean profile
!@var PEDNL00 edge pressure for mean profile (mb)
      REAL*8, DIMENSION(LM+LM_REQ) ::
     &     PL00, PMIDL00, PDSIGL00, AML00, BYAML00
      REAL*8, DIMENSION(LM+LM_REQ+1) :: PEDNL00

!**** Model control parameters:

!@var ij_debug: if i > 0, print out some extra info on bad ij box
      integer, dimension(2) :: ij_debug = (/ 0 , 1 /)

!**** Diagnostic control parameters
!@dbparam Kradia if -1 save data for, if 1|2 do   inst|adj forcing run
      integer :: Kradia=0,iu_rad

!**** Main atmospheric prognostic variables
!@var MA = Air mass per unit area of each layer (kg/m^2)
!@var MAOLD = MA before dynamics used by advection and condensation
!@var U,V east-west, and north-south velocities (m/s)
!@var T potential temperature (referenced to 1 mb) (K)
!@var Q specific humidity (kg water vapor/kg air)
!@var qcl cloud liquid water amount (kg water/kg air)
!@var qci cloud ice water amount (kg water/kg air)
      Real*8,Allocatable,Dimension(:,:,:) :: MA,MAOLD,U,V,T,Q,QCL,QCI

!@var MASUM (kg/m^2) = [column mass per unit area] - MTOP
!@var P surface pressure (hecto-Pascals - PTOP)
      Real*8,Allocatable,Dimension(:,:) :: MASUM,P

      real*8, parameter :: temperature_istart1=250. ! not used

! flag needed until surface components always obtain near-surface conditions
! from the atm state rather than reading them directly from the AIC
      logical :: traditional_coldstart_aic=.true.

!**** Boundary condition arrays:
!@var ZATMO: surface elevation (m)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: ZATMO

C**** Some helpful arrays (arrays should be L first)
!@var  PDSIG  Surface pressure * DSIG(L) (mb)
!@var  byMA = 1/MA (m^2/kg)
!@var  PMID  Pressure at mid point of box (mb)
!@var  PK   PMID**KAPA
!@var  PEDN  Pressure at lower edge of box (incl. surface) (mb)
!@var  PEK  PEDN**KAPA
!@var  PTROPO  Pressure at mid point of tropopause level (mb)
!@var  LTROPO  Tropopause layer
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PDSIG
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: byMA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PMID
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PMIDOLD
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PK
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PEDN
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PEK
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PTROPO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LTROPO
#ifdef etc_subdd
!@var  TTROPO  Temperature at mid point of tropopause level, extra subdaily
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: TTROPO
#endif

C**** module should own dynam variables used by other routines
!@var SD_CLOUDS vert. integrated horizontal convergence (for clouds)
!@var GZ geopotential height (for Clouds and Diagnostics)
!@var DPDX_BY_RHO,DPDY_BY_RHO (pressure gradients)/density at L=1
!@var DPDX_BY_RHO_0,DPDY_BY_RHO_0 surface (pressure gradients)/density
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SD_CLOUDS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GZ
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDX_BY_RHO,DPDY_BY_RHO
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDX_BY_RHO_0
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDY_BY_RHO_0
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PHI

!@var MUs,MVs,MWs,MB save for source time step tracer advection
!@var MB Air mass array for tracers (before advection)
!@var MMA (kg) Air mass array for tracers (updated during advection)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: MUs,MVs,MWs,MB,MMA

!@var DKE change in KE due to dissipation (SURF/DC/MC) (m^2/s^2)
!@var KEA KE on the A grid (m^2/s^2)
!@var UALIJ,VALIJ U,V on the A grid (m/s)
!@var WSAVE vertical velocity (m/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DKE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: KEA ! ke on A grid
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: UALIJ,VALIJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WSAVE

      END MODULE ATM_COM


      SUBROUTINE ALLOC_ATM_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE CONSTANT, only : GRAV
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID,HALO_UPDATE
     &     ,hassouthpole,hasnorthpole
      Use RESOLUTION, Only: IM,JM,LM, MDRYA
      USE ATM_COM, ONLY : temperature_istart1
      USE ATM_COM, ONLY : ZATMO,P,U,V,T,Q,qcl,qci
      USE ATM_COM, ONLY :
     &     PDSIG,MA,MAOLD,byMA,PMID,PMIDOLD,PK,
     &     PEDN,PEK,SD_CLOUDS,GZ,PHI,
     &     MUs,MVs,MWs,MB,MMA,DKE,KEA,
     &     UALIJ,VALIJ,WSAVE,
     &     MASUM,PTROPO,LTROPO,
     &     DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
      use pario, only : par_open,par_close,read_dist_data
      use Dictionary_mod, only : sync_param, get_param
#ifdef etc_subdd
         Use ATM_COM, Only: TTROPO
#endif

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid
      INTEGER :: fid
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: I, I_0, I_1, J_1, J_0
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP

C****
C**** CALCULATE SPHERICAL GEOMETRY
C****
      ALLOCATE(ZATMO(I_0H:I_1H,J_0H:J_1H), STAT = IER)
      ALLOCATE(P(I_0H:I_1H,J_0H:J_1H), STAT = IER)
      ALLOCATE(U(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
      ALLOCATE(V(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
      ALLOCATE(T(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
      ALLOCATE(Q(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
      ALLOCATE(qcl(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
      ALLOCATE(qci(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      U(:,:,:)=0.
      V(:,:,:)=0.
      T(:,:,:)=temperature_istart1  ! will be changed to pot.temp later
      Q(:,:,:)=3.D-6
      !P(:,:)=PSFMPT
      qcl(:,:,:)=0.
      qci(:,:,:)=0.
      ZATMO(:,:)=0.

      fid = par_open(grid,'TOPO','read')
      call read_dist_data(grid,fid,'zatmo',zatmo)
      call par_close(grid,fid)

      ZATMO(I_0:I_1,J_0:J_1) = ZATMO(I_0:I_1,J_0:J_1)*GRAV   ! Geopotential
      CALL HALO_UPDATE(GRID, ZATMO)
C**** Check polar uniformity
      if(hassouthpole(grid)) then
        do i=2,im
          if (zatmo(i,1).ne.zatmo(1,1)) then
            print*,"Polar topography not uniform, corrected",i,1
     *           ,zatmo(i,1),zatmo(1,1)
            zatmo(i,1)=zatmo(1,1)
          end if
        end do
      end if
      if(hasnorthpole(grid)) then
        do i=2,im
          if (zatmo(i,jm).ne.zatmo(1,jm)) then
            print*,"Polar topography not uniform, corrected",i,jm
     *           ,zatmo(i,jm),zatmo(1,jm)
            zatmo(i,jm)=zatmo(1,jm)
          end if
        end do
      end if

!**** Allocate space for (L,I,J) arrays
      ALLOCATE (PDSIG(LM,I_0H:I_1H,J_0H:J_1H),
     $             MA(LM,I_0H:I_1H,J_0H:J_1H),
     $          MAOLD(LM,I_0H:I_1H,J_0H:J_1H),
     $           byMA(LM,I_0H:I_1H,J_0H:J_1H),
     $           PMID(LM,I_0H:I_1H,J_0H:J_1H),
     $        PMIDOLD(LM,I_0H:I_1H,J_0H:J_1H),
     $             PK(LM,I_0H:I_1H,J_0H:J_1H),
     $         PEDN(LM+1,I_0H:I_1H,J_0H:J_1H),
     $          PEK(LM+1,I_0H:I_1H,J_0H:J_1H),
     $   STAT = IER)

      MA(:,:,:) = MDRYA  !  needed for undefined halo cells

!**** Allocate space for (I,J,L) arrays
      ALLOCATE( SD_CLOUDS(I_0H:I_1H,J_0H:J_1H,LM),
     $                 GZ(I_0H:I_1H,J_0H:J_1H,LM),
     $                PHI(I_0H:I_1H,J_0H:J_1H,LM),
     $                MUs(I_0H:I_1H,J_0H:J_1H,LM),
     $                MVs(I_0H:I_1H,J_0H:J_1H,LM),
     $                MWs(I_0H:I_1H,J_0H:J_1H,LM),
     $                 MB(I_0H:I_1H,J_0H:J_1H,LM),
     $                MMA(I_0H:I_1H,J_0H:J_1H,LM),
     $                DKE(I_0H:I_1H,J_0H:J_1H,LM),
     $                KEA(I_0H:I_1H,J_0H:J_1H,LM),
     $                UALIJ(LM,I_0H:I_1H,J_0H:J_1H),
     $                VALIJ(LM,I_0H:I_1H,J_0H:J_1H),
     $              WSAVE(I_0H:I_1H,J_0H:J_1H,LM-1),
     $   STAT = IER)

!**** Allocate space for (I,J) arrays
      Allocate (MASUM(I_0H:I_1H,J_0H:J_1H),
     $          PTROPO(I_0H:I_1H,J_0H:J_1H),
     $          LTROPO(I_0H:I_1H,J_0H:J_1H),
#ifdef etc_subdd
     $          TTROPO(I_0H:I_1H,J_0H:J_1H),   ! extra subdaily
#endif
     $     DPDX_BY_RHO(I_0H:I_1H,J_0H:J_1H),
     $     DPDY_BY_RHO(I_0H:I_1H,J_0H:J_1H),
     $   DPDX_BY_RHO_0(I_0H:I_1H,J_0H:J_1H),
     $   DPDY_BY_RHO_0(I_0H:I_1H,J_0H:J_1H),
     $   STAT = IER)

! correct or wrong, but being static all arrays were initialized
! to zero by default. They have to be initialized to something now
! to avoid floating point exceptions...
      DPDX_BY_RHO(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDY_BY_RHO(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDX_BY_RHO_0(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDY_BY_RHO_0(I_0H:I_1H,J_0H:J_1H) = 0.d0

      SD_CLOUDS(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0

      END SUBROUTINE ALLOC_ATM_COM

ccc was not sure where to dump these routines ... IA
      module conserv_diags
      implicit none

      contains
      subroutine declare_conserv_diags( grid, fid, name_dims )
      use domain_decomp_atm, only : dist_grid, getDomainBounds
      use pario, only : defvar
      implicit none
      type (dist_grid), intent(in) :: grid
      integer ::  fid
      character(len=*) :: name_dims
      integer :: i_0h, i_1h, j_0h, j_1h
      integer :: ier
      real*8, allocatable :: buf(:,:)
      call getDomainBounds( grid, j_strt_halo=j_0h, j_stop_halo=j_1h,
     &     i_strt_halo=i_0h, i_stop_halo=i_1h )
      allocate( buf(i_0h:i_1h,j_0h:j_1h), stat=ier)
      call defvar(grid, fid, buf, trim(name_dims))
      deallocate( buf )
      end subroutine declare_conserv_diags

      subroutine dump_conserv_diags( grid, fid, name, conserv )
      use domain_decomp_atm, only : dist_grid, getDomainBounds
      use pario, only : write_dist_data
      implicit none
      type (dist_grid), intent(in) :: grid
      integer ::  fid
      character(len=*) :: name
      external :: conserv
      integer :: i_0h, i_1h, j_0h, j_1h
      integer :: ier
      real*8, allocatable :: buf(:,:)
      call getDomainBounds( grid, j_strt_halo=j_0h, j_stop_halo=j_1h,
     &     i_strt_halo=i_0h, i_stop_halo=i_1h )
      allocate( buf(i_0h:i_1h,j_0h:j_1h), stat=ier)
      call conserv(buf)
      call write_dist_data(grid,fid,trim(name),buf)
      deallocate( buf )
      end subroutine dump_conserv_diags

      end module conserv_diags


      subroutine def_rsf_atm(fid)
!@sum  def_rsf_model defines U,V,T,P,Q,qcl array structure in restart files
!@auth M. Kelley
!@ver  beta
      use atm_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use conserv_diags
      implicit none
      integer fid   !@var fid file id
      character(len=20) :: ijlstr
      ijlstr='(dist_im,dist_jm,lm)'
      call defvar(grid,fid,u,'u'//ijlstr)
      call defvar(grid,fid,v,'v'//ijlstr)
      call defvar(grid,fid,t,'t'//ijlstr)
      call defvar(grid,fid,q,'q'//ijlstr)
      call defvar(grid,fid,qcl,'qcl'//ijlstr)
      call defvar(grid,fid,qci,'qci'//ijlstr)
      call defvar(grid,fid,p,'p(dist_im,dist_jm)')
      call defvar(grid,fid,ma,'ma(lm,dist_im,dist_jm)')
#ifdef BLK_2MOM
#endif
      call declare_conserv_diags( grid, fid, 'watmo(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'ekatmo(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'epatmo(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'ewatmo(dist_im,dist_jm)' )
      return
      end subroutine def_rsf_atm


      subroutine new_io_atm(fid,iaction)
!@sum  new_io_model read/write U,V,T,P,Q,qcl arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : iowrite,ioread
      use atm_com
      use domain_decomp_atm, only: grid, HALO_UPDATE_COLUMN
      use pario, only : write_dist_data,read_dist_data
      use conserv_diags
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      external conserv_WM, conserv_KE, conserv_PE, conserv_EWM

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 'ma', ma, jdim=3)
        call write_dist_data(grid, fid, 'u', u)
        call write_dist_data(grid, fid, 'v', v)
        call write_dist_data(grid, fid, 't', t)
        call write_dist_data(grid, fid, 'p', p)
        call write_dist_data(grid, fid, 'q', q)
        call write_dist_data(grid, fid, 'qcl', qcl)
        call write_dist_data(grid, fid, 'qci', qci)
#ifdef BLK_2MOM
#endif
        call dump_conserv_diags( grid, fid, 'watmo', conserv_WM )
        call dump_conserv_diags( grid, fid, 'ekatmo', conserv_KE )
        call dump_conserv_diags( grid, fid, 'epatmo', conserv_PE )
        call dump_conserv_diags( grid, fid, 'ewatmo', conserv_EWM  )
      case (ioread)             ! input from restart file
        call read_dist_data(grid, fid, 'ma', ma, jdim=3)
        call read_dist_data(grid, fid, 'u', u)
        call read_dist_data(grid, fid, 'v', v)
        call read_dist_data(grid, fid, 't', t)
        call read_dist_data(grid, fid, 'p', p)
        call read_dist_data(grid, fid, 'q', q)
        call read_dist_data(grid, fid, 'qcl', qcl)
        call read_dist_data(grid, fid, 'qci', qci)
        Call HALO_UPDATE_COLUMN (GRID, MA)
        Call MAtoPMB
#ifdef BLK_2MOM
#endif
      end select
      return
      end subroutine new_io_atm
