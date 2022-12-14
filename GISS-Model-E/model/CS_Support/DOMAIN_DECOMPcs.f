      module domain_decomp_atm
!@auth M. Kelley
!@ver 1.0
!@sum  This module bundles the procedures necessary to a modelE
!@+    cubed-sphere domain decomposition into a single module,
!@+    and provides a version of init_grid specific to a cubed
!@+    sphere atmosphere that runs the GSFC FVcubed dynamical core.

c get grid-independent procedures from domain_decomp_1d
      use domain_decomp_1d, only : am_i_root, sumxpe,
     &     read_parallel,write_parallel,broadcast,
     &     globalmax, setMpiCommunicator,
     &     hasSouthPole, hasNorthPole, getDomainBounds

c get dist_grid, halo_update, globalsum, etc. from the dd2d_utils module
      use dd2d_utils, only : dist_grid,init_dist_grid
     &                      ,halo_update,globalsum
c the following are also available, but use of global arrays is discouraged
c     &                     ,pack_data,unpack_data

c get i/o routines for fortran binary sequential access
      use pario_fbsa  ! readt_parallel et al.

c get i/o routines for netcdf format
      use pario

c get procedures from the MPP package as they become compatible with
c model E or vice versa
c      use mpp_domains_mod

      USE ESMF, only : ESMF_VM,ESMF_Grid

      implicit none
      save

c declare an instance of dist_grid for the atmosphere
      type(dist_grid), target :: grid

c declare requisite instances of these types for the atmosphere
      Type (ESMF_VM) :: modelE_vm
      TYPE (ESMF_Grid) :: ESMF_GRID_ATM

      contains

      SUBROUTINE INIT_GRID(grd_dum, IM,JM,LM,
     &     width,J_SCM,bc_periodic,CREATE_CAP) ! optional arguments
      USE ESMF
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
      INTEGER, INTENT(IN) :: IM,JM,LM
      INTEGER, OPTIONAL :: width
      INTEGER, OPTIONAL, INTENT(IN) :: J_SCM       ! ignored here
      LOGICAL, OPTIONAL, INTENT(IN) :: bc_periodic ! ignored here
      LOGICAL, OPTIONAL, INTENT(IN) :: CREATE_CAP  ! always true for FVcubed
c local variables:
      integer :: rc,width_
      TYPE(ESMF_Grid), external :: AppGridCreateF
      TYPE(ESMF_Config) :: cf
      Type (ESMF_DELayout)::layout
      integer :: npes,npesx,npesy,  rank_glob,rank_tile,rank_i,rank_j
      integer, parameter :: npesx_max=100
      integer, dimension(npesx_max) :: istrt,istop
c      integer :: ntiles, npx, npy, ng, isd, ied , jsd, jed ! for MPP
      include 'mpif.h' ! temporary to see MPI_COMM_WORLD

c sanity checks
      if(jm.ne.im)
     &     call stop_model('cubed sphere init_grid: im != jm',255)
      if(mod(im,2).ne.0)
     &     call stop_model('cubed sphere init_grid: odd im',255)

      width_ = 1
      If (Present(width)) width_=width

      grd_dum%private%PERIODICBC = .false.
      grd_dum%private%hasSouthPole = .false.
      grd_dum%private%hasNorthPole = .false.
      grd_dum%private%hasEquator = .false.

      Call ESMF_Initialize(vm=modelE_vm,
!     &     logkindflag=ESMF_LOGKIND_NONE,
     &     rc=rc)
      Call ESMF_VMGet(modelE_vm, petCount=NPES, localPET=rank_glob)

      npesx = int(floor(sqrt(real(NPES/6))))
      npesy = NPES / npesx  ! in ESMF, npesy = 6*npesx

c FVcubed uses a domain decomposition decided by the MPP library.
c To _temporarily_ avoid build dependencies upon MPP, mimic
c its procedure for grouping IM gridcells into NPESX divisions.
c Store start/end indices for each division in istrt,istop.
      call mimic_mpp_division_of_points(im,npesx,istrt,istop)

      rank_tile = mod(rank_glob,npes/6)
      rank_i    = mod(rank_tile,npesx)
      rank_j    = rank_tile/npesx
      grd_dum%I_STRT        = istrt(1+rank_i)
      grd_dum%I_STOP        = istop(1+rank_i)
      grd_dum%J_STRT        = istrt(1+rank_j)
      grd_dum%J_STOP        = istop(1+rank_j)

c
c Once the MPP communication procedures start to be used, we will
c initialize its domain2D object and retrieve I/J_STRT/STOP from it
c
c      npx=IM+1; npy=JM+1; ng = width_;
c      ntiles=6
c      my_pet = mpp_pe()
c      NPES   = mpp_npes()
c      root   = mpp_root_pe()
c      call init_domain(grd_dum%domain,npx,npy,ntiles,ng,
c     &                 grd_dum%bc_periodic)
c      call mpp_get_compute_domain( grd_dum%domain,
c     &                         I0_DUM, I1_DUM, J0_DUM, J1_DUM)
c      call mpp_get_data_domain( grd_dum%domain, isd, ied , jsd, jed)
c      write(*,*)'mpp-bounds',my_pet,I0_DUM, I1_DUM, J0_DUM, J1_DUM


      grd_dum%IM_WORLD      = IM
      grd_dum%JM_WORLD      = IM

      grd_dum%I_STRT_HALO   = grd_dum%I_STRT - width_
      grd_dum%I_STOP_HALO   = grd_dum%I_STOP + width_
      grd_dum%J_STRT_HALO   = grd_dum%J_STRT - width_
      grd_dum%J_STOP_HALO   = grd_dum%J_STOP + width_

      grd_dum%J_STRT_SKP    = grd_dum%J_STRT 
      grd_dum%J_STOP_SKP    = grd_dum%J_STOP

c set these to something huge to trigger segfaults if accidentlly used?
      grd_dum%J_STRT_STGR   = grd_dum%J_STRT
      grd_dum%J_STOP_STGR   = grd_dum%J_STOP

      call setMpiCommunicator(grd_dum, MPI_COMM_WORLD)

c
c initialize compoments of dist_grid specific to dd2d_utils routines
c
      call init_dist_grid(
     &     grd_dum%IM_WORLD,grd_dum%JM_WORLD,6,
     &     grd_dum%I_STRT,grd_dum%I_STOP,
     &     grd_dum%J_STRT,grd_dum%J_STOP,
     &     grd_dum%I_STRT_HALO,grd_dum%I_STOP_HALO,
     &     grd_dum%J_STRT_HALO,grd_dum%J_STOP_HALO,
     &     MPI_COMM_WORLD, grd_dum)

      grd_dum%have_domain = .true.
      grd_dum%mpi_comm = MPI_COMM_WORLD
      grd_dum%npes_comm = grd_dum%nproc

c
c create the configuration file needed by the FVcubed core
c
      cf = load_cap_config('cap.rc',IM,JM*6,LM,npesx,npesy)
      print*, 'Started AppGridCreateF'
      ESMF_GRID_ATM = AppGridCreateF(cf, modelE_vm, rc)
c is this needed?
      !call ESMF_GridGet(ESMF_GRID_ATM, delayout=layout, rc=rc)
      print*, 'Finished AppGridCreateF'

      return
      END SUBROUTINE INIT_GRID

      subroutine mimic_mpp_division_of_points(npx,ndiv,is,ie)
c Divides 1:npx gridpoints into ndiv intervals.
c The start,end indices of interval m are stored in is(m),ie(m).
c Intended only for even values of npx!
      implicit none
      integer :: npx,ndiv ! inputs
      integer, dimension(ndiv) :: is,ie ! outputs
      integer :: idiv,ni,jdiv
      real*4 :: rhalf
      rhalf = .5*real(ndiv)
      is(1) = 1
      ie(ndiv) = npx
      do idiv=1,ndiv/2      ! loop over half of the intervals
        jdiv = ndiv-idiv+1  ! and set the other half using symmetry
        if(idiv.gt.1) then
          is(idiv) = ie(idiv-1) + 1
          ie(jdiv) = is(jdiv+1) - 1
        endif
! the following formula for ni is based on the one in MPP:
! domain is sized by dividing the remaining points by remaining domains
! CEILING causes the larger intervals to lie at the edges
        ni = CEILING( real(npx/2-is(idiv)+1)/(rhalf-real(idiv-1)) )
        ie(idiv) = is(idiv) + ni - 1
        is(jdiv) = ie(jdiv) - ni + 1
      enddo
      if(ndiv.gt.1 .and. mod(ndiv,2).eq.1) then
c ndiv odd: the middle interval
        idiv = 1+ndiv/2
        is(idiv) = ie(idiv-1) + 1
        ie(idiv) = is(idiv+1) - 1
      endif
      return
      end subroutine mimic_mpp_division_of_points

! ----------------------------------------------------------------------
      function load_cap_config(config_file,IM,JM,LM,NP_X,NP_Y)
     &     result( config )
! ----------------------------------------------------------------------
      use ESMF
      use FILEMANAGER    
      include 'mpif.h' ! temporary to see MPI_Barrier
      character(len=*), parameter :: Iam=
     &     "DOMAIN_DECOMP::load_cap_config"
      character(len=*), intent(in) :: config_file
      integer,          intent(in) :: IM,JM,LM,NP_X,NP_Y
      type (esmf_config)           :: config
    
      integer :: rc, iunit
    
      config = ESMF_ConfigCreate(rc=rc)
      call openunit(config_file, iunit, qbin=.false., qold=.false.)
      if(am_i_root()) then
        write(iunit,*)'IM:  ', IM
        write(iunit,*)'JM:  ', JM
        write(iunit,*)'LM:  ', LM
        write(iunit,*)'NX:  ', NP_X
        write(iunit,*)'NY:  ', NP_Y
      endif
      call closeUnit(iunit)
     
      call MPI_Barrier(mpi_comm_world, rc)
      call ESMF_ConfigLoadFile(config, config_file, rc=rc)

      end function load_cap_config

      end module domain_decomp_atm
