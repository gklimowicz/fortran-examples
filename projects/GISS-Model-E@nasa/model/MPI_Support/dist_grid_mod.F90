#define VERIFY_(rc) If (rc /= ESMF_SUCCESS) Call abort(__LINE__,rc)
#include "rundeck_opts.h"
#ifdef MPI_DEFS_HACK
#include "mpi_defs.h"
#endif

#ifndef CUBED_SPHERE
#define DOMAIN_DECOMP_ATM_IS_1D
#endif

#ifdef NEW_IO
#define USE_DD2D_UTILS
#endif

#ifdef CUBED_SPHERE
#define USE_DD2D_UTILS
#undef CUBED_SPHERE
#undef USE_ESMF /* this is handled elsewhere for the moment */
#endif

MODULE dist_grid_mod
!@sum  DOMAIN_DECOMP encapsulates lat-lon decomposition information
!@+    for the message passing (ESMF) implementation.
!@auth NCCS ASTG

! NOTE: USE_ESMF and USE_MPI are #defined from the "make" command line
! using options ESMF=YES and MPI=YES respectively. When ESMF=YES, MPI=YES
! is also implied. However MPI=YES excludes ESMF.

#ifdef USE_ESMF
   use ESMF
#endif
   use MpiSupport_mod, only: am_i_root
   use Domain_mod
   use Hidden_mod

#ifdef USE_DD2D_UTILS
   use dd2d_utils, only : dist_grid,init_dist_grid
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif
   SAVE
   PRIVATE ! Except for

!@var DIST_GRID derived type to provide domain decomposition information
!@+   Public components are used to minimize overhead for accessing
!@+   routine components
   PUBLIC :: DIST_GRID
   public :: setMpiCommunicator
   public :: setCommunicator
!@var INIT_APP set some parameters and initialize ESMF
   PUBLIC :: INIT_APP
   PUBLIC :: INIT_GRID
   PUBLIC :: DESTROY_GRID
!@var FINISH_APP Cleans up at the end of the run (closes debugging file)
   PUBLIC :: FINISH_APP
!@var GLOBALMIN determine max value across pes
   PUBLIC :: GLOBALMIN
!@var GLOBALMAX determine max value across pes
   PUBLIC :: GLOBALMAX
!@var SUMXPE sum an array over processors without reducing its rank
   PUBLIC :: SUMXPE
!@var GET - extracts bounds information from DIST_GRID object
   PUBLIC :: getDomainBounds
   PUBLIC :: HERE
   PUBLIC :: LOG_PARALLEL
  
   PUBLIC :: TRANSP
   PUBLIC :: TRANSPOSE_COLUMN
!@var GLOBALMIN Generic wrapper for Real
   INTERFACE GLOBALMIN
      MODULE PROCEDURE GLOBALMIN_I
      MODULE PROCEDURE GLOBALMIN_R
   END INTERFACE

!@var GLOBALMAX Generic wrapper for Real/integer
   INTERFACE GLOBALMAX
      MODULE PROCEDURE GLOBALMAX_R
      MODULE PROCEDURE GLOBALMAX_I
      MODULE PROCEDURE GLOBALMAX_I_1D
   END INTERFACE

   INTERFACE TRANSP
      MODULE PROCEDURE TRANSPOSE_ij
      MODULE PROCEDURE TRANSPOSE_ijk
   END INTERFACE
  
   INTERFACE SUMXPE
      MODULE PROCEDURE SUMXPE_1D
      MODULE PROCEDURE SUMXPE_1D_I
      MODULE PROCEDURE SUMXPE_2D
      MODULE PROCEDURE SUMXPE_3D
      MODULE PROCEDURE SUMXPE_4D
!        MODULE PROCEDURE SUMXPE_5D
   END INTERFACE

!@var broadcast Generic routine to broadcast data to all PEs.
   PUBLIC :: broadcast
   INTERFACE broadcast
      MODULE PROCEDURE broadcast_0D
      MODULE PROCEDURE broadcast_1D
      MODULE PROCEDURE broadcast_2D
      MODULE PROCEDURE broadcast_3D
      MODULE PROCEDURE broadcast_4D
      MODULE PROCEDURE ibroadcast_0D
      MODULE PROCEDURE ibroadcast_0D_world
      MODULE PROCEDURE ibroadcast_1D
      MODULE PROCEDURE ibroadcast_2D
      MODULE PROCEDURE ibroadcast_3D
      MODULE PROCEDURE ibroadcast_4D
   END INTERFACE

!@var BAND_PACK Procedure in which each PE receives data from other PEs
!@+             to fill a contiguous pre-requested range of J indices
   public :: band_pack,band_pack_column
   interface band_pack
      module procedure band_pack_ij
      module procedure band_pack_ijl
   end interface
!@var BAND_PACK_TYPE a data structure needed by BAND_PACK, initialized
!@var via INIT_BAND_PACK_TYPE
   type band_pack_type
      integer :: j_strt,j_stop
      integer :: j_strt_halo,j_stop_halo
      integer :: jband_strt,jband_stop
      integer, dimension(:), pointer :: scnts,sdspl,sdspl_inplace
      integer, dimension(:), pointer :: rcnts,rdspl,rdspl_inplace
      integer, dimension(:), pointer :: j0_send,j1_send
      integer, dimension(:), pointer :: j0_recv,j1_recv
      integer :: mpi_comm
      integer :: npes_comm
   end type band_pack_type
!@var INIT_BAND_PACK_TYPE initialization routine during which each PE
!@+   requests a range of J indices and sets up the necessary send/receive
!@+   information for the BAND_PACK procedure
   public :: band_pack_type,init_band_pack_type

   PUBLIC SEND_TO_J
   interface SEND_TO_J
      module procedure SEND_TO_J_1D
      module procedure ISEND_TO_J_0D
   end interface

   PUBLIC RECV_FROM_J
   interface RECV_FROM_J
      module procedure RECV_FROM_J_1D
      module procedure IRECV_FROM_J_0D
   end interface

   INTEGER, PARAMETER :: HALO_WIDTH = 1
   integer ::  root

#ifndef USE_DD2D_UTILS
      ! Local grid information
   TYPE DIST_GRID

      type (Hidden_type) :: private
#ifdef USE_ESMF
      TYPE (ESMF_Grid) :: ESMF_GRID
#endif
      INTEGER :: NPES_USED
      ! Parameters for Global domain
      INTEGER :: IM_WORLD        ! Number of Longitudes
      INTEGER :: JM_WORLD        ! Number of latitudes
      ! Parameters for local domain
      LOGICAL :: HAVE_DOMAIN = .false.     ! Whether this PE has any of the domain
      INTEGER :: I_STRT          ! Begin local domain longitude index
      INTEGER :: I_STOP          ! End   local domain longitude index
      INTEGER :: J_STRT          ! Begin local domain latitude  index
      INTEGER :: J_STOP          ! End   local domain latitude  index
      INTEGER :: J_STRT_SKP      ! Begin local domain exclusive of S pole
      INTEGER :: J_STOP_SKP      ! End   local domain exclusive of N pole
      INTEGER :: ni_loc ! for transpose
      ! Parameters for halo of local domain
      INTEGER :: I_STRT_HALO     ! Begin halo longitude index
      INTEGER :: I_STOP_HALO     ! End   halo longitude index
      INTEGER :: J_STRT_HALO     ! Begin halo latitude  index
      INTEGER :: J_STOP_HALO     ! End   halo latitude  index
      ! Parameters for staggered "B" grid
      ! Note that global staggered grid begins at "2".
      INTEGER :: J_STRT_STGR     ! Begin local staggered domain
      INTEGER :: J_STOP_STGR     ! End   local staggered domain
      
      INTEGER, DIMENSION(:), POINTER :: DJ_MAP
      INTEGER :: DJ
      
   END TYPE DIST_GRID
#endif

   public :: haveLatitude
   public :: isinLocalSubdomain

   public :: rank, npes_world
   INTEGER :: communicator
!@var NPES_WORLD number of total processes
   INTEGER :: NPES_WORLD
!@var RANK index of _this_ PET (analagous to MPI rank)
   INTEGER :: rank

   INTEGER, PUBLIC :: CHECKSUM_UNIT

   public :: barrier
   public :: isPeriodic
   public :: am_i_root
   ! Direction bits
   public :: NORTH, SOUTH, NORTH2, NORTH3, SOUTHJMm1
   Integer, Parameter :: NORTH  = 2**0
   Integer, Parameter :: NORTH2 = 2**2
   Integer, Parameter :: NORTH3 = 2**3
   Integer, Parameter :: SOUTH  = 2**1
   Integer, Parameter :: SOUTHJMm1  = 2**4
   public :: getMpiCommunicator
   public :: getNumProcesses
   public :: getNumAllProcesses
   public :: getMpiTag
   public :: incrementMpiTag
   public :: hasSouthPole
   public :: hasNorthPole
   public :: hasPeriodicBC  
   public :: getLogUnit

#ifdef USE_ESMF
   public :: load_cap_config
   Public :: modelE_vm
   Type (ESMF_VM), Target :: modelE_vm
#endif

#ifdef USE_MPI
   type axisIndex
      sequence  ! sequence forces the data elements 
                ! to be next to each other in memory 
      integer :: min
      integer :: max
   end type axisIndex
   Public :: axisIndex 
   Public :: getAxisIndex
   interface getAxisIndex
#ifdef USE_ESMF
      module procedure getESMFAxisIndex
#endif
#ifndef USE_ESMF
      module procedure getMPIAxisIndex
#endif
   end interface
   public :: gridRootPELocation
   public :: gridPELayout
#endif

   integer, parameter :: maxStrLen = 40
   public :: maxStrLen

 CONTAINS

! ----------------------------------------------------------------------
   ! This routine initializes the quantities described above.
   ! The initialization should proceed prior to any grid computations.
   subroutine init_app
! ----------------------------------------------------------------------
     USE FILEMANAGER, ONLY : openunit
     integer :: rc, comm

     rank = 0        ! default rank = root PE for serial run
     NPES_WORLD = 1    ! default NPES = 1 for serial run
#ifdef USE_ESMF
   call ESMF_Initialize(vm=modelE_vm, logkindflag=ESMF_LOGKIND_NONE, rc=rc)
   VERIFY_(rc)
   call ESMF_VMGet(modelE_vm, localPET=rank, petCount=NPES_WORLD, &
    &     mpiCommunicator=comm, rc=rc)
   VERIFY_(rc)
   call setCommunicator(comm)
#else
#ifdef USE_MPI
   call MPI_INIT(rc)
   call setCommunicator(MPI_COMM_WORLD)
   call MPI_COMM_SIZE(COMMUNICATOR, NPES_WORLD, rc)
   call MPI_COMM_RANK(COMMUNICATOR, rank, rc)
#endif
#endif
     if (rank == 0) write(*,*)'Num MPI Processes: ', NPES_WORLD
#ifdef DEBUG_DECOMP
     CHECKSUM_UNIT = rank
     CALL openunit('debug_decomp',CHECKSUM_UNIT)
#endif
     return
   end subroutine init_app


! ----------------------------------------------------------------------
   subroutine INIT_GRID(distGrid, IM, JM, LM, &
        width, J_SCM, bc_periodic, CREATE_CAP, npes_max)
! ----------------------------------------------------------------------
     USE FILEMANAGER, Only : openunit
     IMPLICIT NONE
     TYPE (DIST_GRID), INTENT(INOUT) :: distGrid
     INTEGER, INTENT(IN) :: IM, JM,LM
     INTEGER, OPTIONAL, INTENT(IN) :: J_SCM ! single column model
     INTEGER, OPTIONAL :: width
     LOGICAL, OPTIONAL, INTENT(IN) :: bc_periodic
     LOGICAL, OPTIONAL, INTENT(IN) :: CREATE_CAP
     INTEGER, OPTIONAL, INTENT(IN) :: npes_max
#ifdef USE_MPI
     Type (AxisIndex), Pointer :: AI(:,:)
     integer, allocatable   :: IMS(:), JMS(:)
#endif
     INTEGER :: p, rc, ierr
     integer :: npes_used
     integer, dimension(:), allocatable :: pelist
     INTEGER :: J_EQUATOR
     INTEGER :: width_
     INTEGER :: AIbounds(4)
     integer :: group_world, group_used
     integer :: newCommunicator   
     character(len=256) :: errorMessage  

     distGrid%IM_WORLD      = IM
     distGrid%JM_WORLD      = JM

#ifdef USE_MPI
     call MPI_COMM_SIZE(COMMUNICATOR, distGrid%NPES_WORLD, rc)
     call MPI_COMM_RANK(COMMUNICATOR, distGrid%rank, rc)
! repeated assignment is needed for MPI unit tests
     NPES_WORLD = distGrid%NPES_WORLD
     rank = distGrid%rank
     Allocate(distGrid%dj_map(0:npes_world-1))

    ! Distribute the horizontal dimensions
    !-------------------------------------
     allocate(ims(0:0), jms(0:distGrid%npes_world-1))
     distGrid%npes_used = min(distGrid%npes_world, jm)
     if (present(npes_max))distGrid% npes_used = min(npes_max, distGrid%npes_used)
     jms(0:distGrid%npes_used-1) = getLatitudeDistribution(distGrid)
     if(distGrid%npes_used<distGrid%npes_world) &
       jms(distGrid%npes_used:distGrid%npes_world-1) = 0
     ims(0) = im
#ifndef USE_ESMF
     AIbounds = MPIgridBounds(distGrid)
#endif

#else  /* SERIAL CASE */

! repeated assignment is needed for MPI unit tests
     rank = 0        ! default rank = root PE for serial run
     NPES_WORLD = 1    ! default NPES = 1 for serial run
     distGrid%rank = 0
     distGrid%npes_world = 1
     distGrid%npes_used = 1
     AIbounds = (/ 1, IM, 1, JM /)
     if (present(J_SCM)) then
        AIbounds(3) = J_SCM
        AIbounds(4) = J_SCM
     end if

#endif /* USE_MPI */

#ifdef USE_ESMF
     distGrid%ESMF_GRID = ESMF_GridCreate( &
            name="modelE grid",            &
            countsPerDEDim1=ims,           &
            countsPerDEDim2=jms,           &
            indexFlag = ESMF_INDEX_USER,   &
            gridMemLBound = (/1,1/),       &
            gridEdgeLWidth = (/0,0/),      &
            gridEdgeUWidth = (/0,0/),      &
            coordDep1 = (/1,2/),           &
            coordDep2 = (/1,2/),           &
            rc=rc)
    VERIFY_(rc)
    AIbounds =  ESMFgridBounds(distGrid)
#endif

     width_ = HALO_WIDTH
     If (Present(width)) width_=width

     ! Wrapped ESMF grid
     distGrid%I_STRT        = AIbounds(1)
     distGrid%I_STOP        = AIbounds(2)
     distGrid%I_STRT_HALO   = MAX( 1, AIbounds(1)-width_)
     distGrid%I_STOP_HALO   = MIN(IM, AIbounds(2)+width_)
     distGrid%ni_loc = (distGrid%RANK+1)*IM/distGrid%NPES_used - &
                        distGrid%RANK   *IM/distGrid%NPES_used
     
     distGrid%j_strt        = AIbounds(3)
     distGrid%J_STOP        = AIbounds(4)

     call barrier()

     distGrid%HAVE_DOMAIN   = AIbounds(3) <= JM

     distGrid%J_STRT_SKP = max (   2, AIbounds(3))
     distGrid%J_STOP_SKP = min (JM-1, AIbounds(4))

#ifdef USE_MPI
     distGrid%J_STRT_HALO   = AIbounds(3) - width_
     distGrid%J_STOP_HALO   = AIbounds(4) + width_
     distGrid%private%numProcesses = distGrid%npes_used
     distGrid%private%numAllProcesses = distGrid%npes_world
     distGrid%private%mpi_tag = 10  ! initial value

! Create a new MPI communicator including all the PEs with a nonzero
! domain size.  Even when NPES_USED == NPES_WORLD, this is convenient for
! avoiding collisions of MPI tag sequences.

     call mpi_comm_group(COMMUNICATOR,group_world,ierr)
     allocate(pelist(0:distGrid%npes_used-1))
     do p=0,distGrid%npes_used-1
        pelist(p) = p
     enddo
     call mpi_group_incl(group_world,distGrid%npes_used,pelist,group_used,ierr)
     deallocate(pelist)
     call mpi_comm_create(COMMUNICATOR,group_used, newCommunicator, ierr)
     if(.not. distGrid%HAVE_DOMAIN) newCommunicator = MPI_COMM_NULL
     call setMpiCommunicator(distGrid, newCommunicator)
#else
     distGrid%J_STRT_HALO = MAX(1,  distGrid % J_STRT)
     distGrid%J_STOP_HALO = MIN(JM, distGrid % J_STOP)
#endif
     
     distGrid%J_STRT_STGR   = max(2,AIbounds(3))
     distGrid%J_STOP_STGR   = AIbounds(4)

     distGrid%private%hasSouthPole = AIbounds(3) == 1 
     distGrid%private%hasNorthPole = AIbounds(4) == JM .and. AIbounds(3) <= JM

     J_EQUATOR = JM/2
     distGrid%private%hasEquator =  &
          ( (AIbounds(3) <= J_EQUATOR) .AND. (AIbounds(4) >= J_EQUATOR))

     if (JM==1) then
        distGrid%private%hasSouthPole = .false.
        distGrid%private%hasNorthPole = .false.
        distGrid%private%hasEquator    = .false.
        distGrid%J_STRT_SKP = 1
        distGrid%J_STOP_SKP = 1
     endif

#ifdef USE_DD2D_UTILS
! need to initialize the dd2d version of dist_grid for I/O
     call init_dist_grid( &
    &     distGrid%IM_WORLD,distGrid%JM_WORLD,1,  &
    &     distGrid%I_STRT,distGrid%I_STOP, &
    &     distGrid%j_strt,distGrid%J_STOP, &
    &     distGrid%I_STRT_HALO,distGrid%I_STOP_HALO, &
    &     distGrid%J_STRT_HALO,distGrid%J_STOP_HALO, &
    &     COMMUNICATOR, &
    &     distGrid)
#endif

     if (present(J_SCM)) then
        ! assume J_SCM is in "general position"
        distGrid%private%hasSouthPole = .false.
        distGrid%private%hasNorthPole = .false.
        distGrid%private%hasEquator    = .false.
     endif

     distGrid % private%PERIODICBC = isPeriodic(bc_periodic)

     ! assumption: decomposition along "east-west" direction
     ! not used for lat-lon grids
     if(distGrid % private%PERIODICBC) then
        distGrid%private%hasSouthPole = .false.
        distGrid%private%hasNorthPole = .false.
        distGrid%private%hasEquator    = .false.
     endif

     ! set lookup table PET(J)
     Allocate(distGrid%private%lookup_pet(1:JM))
     distGrid%private%lookup_pet(:) = 0

#ifdef USE_MPI
     ALLOCATE(AI(0:distGrid%npes_world-1,3))
     call getAxisIndex(distGrid, AI)

     Do p = 1, distGrid%npes_world
       distGrid%private%lookup_pet( AI(p-1,2)%min : AI(p-1,2)%max ) = p-1
     end do
     Do p = 0, distGrid%npes_world-1
       distGrid%dj_map(p) = AI(p,2)%max - AI(p,2)%min + 1
     end do
     distGrid%dj=distGrid%dj_map(distGrid%rank)
     
     deallocate(AI)
#endif

   END SUBROUTINE INIT_GRID

! ----------------------------------------------------------------------
   SUBROUTINE DESTROY_GRID(distGrid)
! ----------------------------------------------------------------------
     TYPE (DIST_GRID), INTENT(INOUT) :: distGrid
     integer :: ier
#ifdef USE_ESMF
     Call ESMF_GridDestroy(distGrid%ESMF_Grid)
#endif
#ifdef USE_MPI
     call MPI_Comm_free(distGrid%private%MPI_COMM, ier)
#endif
   END SUBROUTINE DESTROY_GRID

! ----------------------------------------------------------------------
   SUBROUTINE getDomainBounds(distGrid, I_STRT, I_STOP, &
  &                        I_STRT_HALO, I_STOP_HALO, &
  &                        J_STRT, J_STOP, J_STRT_HALO, J_STOP_HALO, &
  &                        J_STRT_SKP, J_STOP_SKP, &
  &                        J_STRT_STGR, J_STOP_STGR, &
  &                        have_south_pole, have_north_pole)
! ----------------------------------------------------------------------
     TYPE (DIST_GRID), INTENT(IN) :: distGrid
     INTEGER, OPTIONAL :: I_STRT, I_STOP
     INTEGER, OPTIONAL :: I_STRT_HALO, I_STOP_HALO
     INTEGER, OPTIONAL :: J_STRT, J_STOP
     INTEGER, OPTIONAL :: J_STRT_HALO, J_STOP_HALO
     INTEGER, OPTIONAL :: J_STRT_SKP, J_STOP_SKP
     INTEGER, OPTIONAL :: J_STRT_STGR, J_STOP_STGR
     LOGICAL, OPTIONAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

     IF (PRESENT(I_STRT)) I_STRT = distGrid%I_STRT
     IF (PRESENT(I_STOP)) I_STOP = distGrid%I_STOP
     
     IF (PRESENT(I_STRT_HALO)) I_STRT_HALO = distGrid%I_STRT_HALO
     IF (PRESENT(I_STOP_HALO)) I_STOP_HALO = distGrid%I_STOP_HALO
     
     IF (PRESENT(J_STRT)) J_STRT = distGrid%j_strt
     IF (PRESENT(J_STOP)) J_STOP = distGrid%J_STOP

     IF (PRESENT(J_STRT_HALO)) J_STRT_HALO = distGrid%J_STRT_HALO
     IF (PRESENT(J_STOP_HALO)) J_STOP_HALO = distGrid%J_STOP_HALO
     
     IF (PRESENT(J_STRT_SKP)) J_STRT_SKP = distGrid%J_STRT_SKP
     IF (PRESENT(J_STOP_SKP)) J_STOP_SKP = distGrid%J_STOP_SKP
     
     IF (PRESENT(J_STRT_STGR)) J_STRT_STGR = distGrid%J_STRT_STGR
     IF (PRESENT(J_STOP_STGR)) J_STOP_STGR = distGrid%J_STOP_STGR
     
     IF (PRESENT(HAVE_SOUTH_POLE)) &
    &             HAVE_SOUTH_POLE= distGrid%private%hasSouthPole
     IF (PRESENT(HAVE_NORTH_POLE)) &
    &             HAVE_NORTH_POLE= distGrid%private%hasNorthPole

   END SUBROUTINE getDomainBounds


! ----------------------------------------------------------------------
   SUBROUTINE FINISH_APP()
! ----------------------------------------------------------------------
     USE FILEMANAGER, ONLY : closeunit
     IMPLICIT NONE
     
     INTEGER :: rc

#ifdef DEBUG_DECOMP
     CALL closeunit(CHECKSUM_UNIT)
     CALL closeunit(grid%private%log_unit)
#endif

#ifdef USE_ESMF
     CALL ESMF_FINALIZE(rc=rc)
#else
#ifdef USE_MPI
   call MPI_Finalize(rc)
#endif
#endif

   END SUBROUTINE FINISH_APP

#ifdef USE_ESMF

! ----------------------------------------------------------------------
   function ESMFgridBounds(grid) result(AIbounds)
! ----------------------------------------------------------------------
     type (DIST_Grid), intent(in) :: grid
     integer :: AIbounds(4)
     type(AxisIndex), dimension(:,:), pointer :: AI

     allocate(AI(0:grid%NPES_WORLD-1,3))
     call getESMFAxisIndex(grid, AI)

     AIbounds(1) = AI(grid%rank,1)%min
     AIbounds(2) = AI(grid%rank,1)%max
     AIbounds(3) = AI(grid%rank,2)%min
     AIbounds(4) = AI(grid%rank,2)%max

     deallocate(AI)

   end function ESMFgridBounds

! ----------------------------------------------------------------------
  subroutine getESMFAxisIndex(grid, AI)
! ----------------------------------------------------------------------
     type (DIST_Grid), intent(in) :: grid
     type (AxisIndex) :: AI(0:,:)
     ! local vars
     integer                               :: status
     character(len=maxStrLen)              :: IAm='getAxisIndex'     
     type (ESMF_DistGrid)                  :: distGrid
     type(ESMF_DELayout)                   :: LAYOUT
     integer,               allocatable    :: AL(:,:)
     integer,               allocatable    :: AU(:,:)
     integer                               :: nDEs
     integer                               :: gridRank
     integer :: p, npes_end, deId, I1,IN,J1,JN
     integer                               :: deList(1)
     
     call ESMF_GridGet(GRID%ESMF_Grid, dimCount=gridRank, distGrid=distGrid, &
          rc=STATUS)
     VERIFY_(status)
     call ESMF_DistGridGet(distGRID, delayout=layout, &
          rc=STATUS)
     VERIFY_(status)
     call ESMF_DELayoutGet(layout, deCount =nDEs, localDeList=deList, &
          rc=status)
     VERIFY_(status)
     deId = deList(1)
     
     allocate (AL(gridRank,0:nDEs-1),  stat=status)
     allocate (AU(gridRank,0:nDEs-1),  stat=status)
     
     call ESMF_DistGridGet(distgrid, minIndexPDe=AL, maxIndexPDe=AU, &
          rc=status)
     VERIFY_(status)
    
     I1 = AL(1, deId)
     IN = AU(1, deId)
     J1 = AL(2, deId)
     JN = AU(2, deId)
     
     do p = 0, nDEs - 1
        AI(p,1)%min = AL(1, p) ! I1
        AI(p,1)%max = AU(1, p) ! IN
        AI(p,2)%min = AL(2, p) ! J1
        AI(p,2)%max = AU(2, p) ! JN
     end do
     
     npes_end = size(AI,1)
     AI(nDEs:npes_end-1,2)%min = AI(nDEs-1,2)%max + 1
     AI(nDEs:npes_end-1,2)%max = AI(nDEs-1,2)%max
     
     deallocate(AU, AL)
     
   end subroutine getESMFAxisIndex

! ----------------------------------------------------------------------
   function load_cap_config(config_file,IM,JM,LM,NP_X,NP_Y)  &
  &         result( config )
! ----------------------------------------------------------------------
     use ESMF
     use FILEMANAGER    
     character(len=*), parameter :: Iam= &
    &               "DOMAIN_DECOMP::load_cap_config"
     character(len=*), intent(in) :: config_file
     integer,          intent(in) :: IM,JM,LM,NP_X,NP_Y
     type (esmf_config)           :: config
    
     integer :: rc, iunit
    
     config = ESMF_ConfigCreate(rc=rc)
     VERIFY_(rc)
     call openunit(config_file, iunit, qbin=.false., qold=.false.)
     if(am_i_root()) then
        write(iunit,*)'IM:  ', IM
        write(iunit,*)'JM:  ', JM
        write(iunit,*)'LM:  ', LM
        write(iunit,*)'NX:  ', NP_X
        write(iunit,*)'NY:  ', NP_Y
     endif
     call closeUnit(iunit)
     
     call barrier()
     call ESMF_ConfigLoadFile(config, config_file, rc=rc)
     
   end function load_cap_config

#else

#ifdef USE_MPI

! ----------------------------------------------------------------------
     function getLatitudeDistribution(grid) result(latsPerProcess)
! ----------------------------------------------------------------------
     TYPE (DIST_GRID), INTENT(IN) :: grid
! Constraint: assumes grid%JM_WORLD >=4.
       integer :: latsPerProcess(0:grid%npes_used-1)
       
       integer :: excess, npes_used, p
       integer :: localAdjustment

       latsPerProcess = 0
       npes_used = grid%npes_used

       ! Set minimum requirements per processor
       ! Currently this is 1 lat/proc away from poles
       ! and 2 lat/proc at poles
       select case (grid%npes_world)
       case (1)
          latsPerProcess = grid%JM_WORLD
          return
       case (2)
          latsPerProcess(0) = grid%JM_WORLD/2
          latsPerProcess(1) = grid%JM_WORLD - (grid%JM_WORLD/2)
          return
       case (3:)
          ! 1st cut - round down
          latsPerProcess(0:npes_used-1) = GRID%JM_WORLD/npes_used 
          
          ! Fix at poles
          ! redistrute excess
          excess = GRID%JM_WORLD - sum(latsPerProcess(0:npes_used-1))
          if (excess > 0) then
            latsPerProcess(0) = max(2, latsPerProcess(0))
          end if
          if (excess > 1) then
            latsPerProcess(npes_used-1) = max(2, latsPerProcess(npes_used-1))
          end if          
          ! re-calculate excess again, since latsPerProcess may have changed
          excess = GRID%JM_WORLD - sum(latsPerProcess(0:npes_used-1))
          ! redistribute any remaining excess among interior processors
          do p = 1, npes_used - 2
             localAdjustment = (p+1)*excess/(npes_used-2) - (p*excess)/(npes_used-2)
             latsPerProcess(p) = latsPerProcess(p) + localAdjustment
          end do
       end select
       
     end function getLatitudeDistribution

! ----------------------------------------------------------------------
   function MPIgridBounds(grid) result(AIbounds)
! ----------------------------------------------------------------------
     type (DIST_Grid), intent(in) :: grid
     integer :: AIbounds(4)
     type(AxisIndex), dimension(:,:), pointer :: AI

     allocate(AI(0:grid%npes_world-1,3))
     call getMPIAxisIndex(grid, AI)

     AIbounds(1) = AI(grid%rank,1)%min
     AIbounds(2) = AI(grid%rank,1)%max
     AIbounds(3) = AI(grid%rank,2)%min
     AIbounds(4) = AI(grid%rank,2)%max

     deallocate(AI)

   end function MPIgridBounds

! ----------------------------------------------------------------------
   subroutine getMPIAxisIndex(grid, AI)
! ----------------------------------------------------------------------
     type (DIST_Grid), intent(in) :: grid
     type (axisIndex) :: AI(0:,:)

     character(len=maxStrLen) :: IAm='getAxisIndex'     
     integer :: p, npes_end
     integer, allocatable   :: jms(:)

     allocate(jms(0:grid%npes_world-1))

     jms(0:grid%npes_used-1) = getLatitudeDistribution(grid)
     if(grid%npes_used<grid%npes_world) jms(grid%npes_used:grid%npes_world-1) = 0

     AI = computeAxisIndex(grid, jms)

     npes_end = size(AI,1)
     AI(grid%npes_world:npes_end-1,2)%min = AI(grid%npes_world-1,2)%max + 1
     AI(grid%npes_world:npes_end-1,2)%max = AI(grid%npes_world-1,2)%max

     deallocate(jms)

   end subroutine getMPIAxisIndex

! ----------------------------------------------------------------------
   function computeAxisIndex(grid, jmsGlob) result(thisAI)
! ----------------------------------------------------------------------
     type (DIST_Grid), intent(in) :: grid
     integer, intent(in) :: jmsGlob(0:)
     type (axisIndex) :: thisAI(0:grid%npes_world-1,3)
     integer :: p 
     
     do p = 0, grid%npes_world - 1
        thisAI(p,1)%min = 1
        thisAI(p,1)%max = grid%IM_WORLD
        if (p==0) then
           thisAI(p,2)%min = 1
           thisAI(p,2)%max = jmsGlob(p)
        else
           thisAI(p,2)%min = thisAI(p-1,2)%max + 1
           thisAI(p,2)%max = thisAI(p-1,2)%max + jmsGlob(p)
        end if
     end do
     
   end function computeAxisIndex

#endif


#endif
      
! ----------------------------------------------------------------------
      SUBROUTINE SUMXPE_1D(arr, arr_master, increment)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:) :: arr
      REAL*8, DIMENSION(:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_MPI
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root(COMMUNICATOR)) then
           allocate(arr_tmp(arr_size))
        else
           allocate(arr_tmp(1))
        end if
        call MPI_Reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,COMMUNICATOR, ierr)
        if(am_i_root(COMMUNICATOR)) then
          arr_master = arr_master + arr_tmp
        endif
        deallocate(arr_tmp)
      else
        call MPI_Reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,COMMUNICATOR, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr is overwritten by itself after reduction
#ifdef USE_MPI
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call MPI_Reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        COMMUNICATOR, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_1D

! ----------------------------------------------------------------------
      SUBROUTINE SUMXPE_1D_I(arr, arr_master, increment)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, DIMENSION(:) :: arr
      INTEGER, DIMENSION(:), optional :: arr_master
      logical, intent(in), optional :: increment
      INTEGER, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_MPI
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root(COMMUNICATOR)) allocate(arr_tmp(arr_size))
        call MPI_Reduce(arr,arr_tmp,arr_size,MPI_INTEGER, &
    &       MPI_SUM,root,COMMUNICATOR, ierr)
        if(am_i_root(COMMUNICATOR)) then
          arr_master = arr_master + arr_tmp
          deallocate(arr_tmp)
        endif
      else
        call MPI_Reduce(arr,arr_master,arr_size,MPI_INTEGER, &
    &       MPI_SUM,root,COMMUNICATOR, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr is overwritten by itself after reduction
#ifdef USE_MPI
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call MPI_Reduce(arr,arr_tmp,arr_size, &
    &        MPI_INTEGER,MPI_SUM,root, &
    &        COMMUNICATOR, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_1D_I

! ----------------------------------------------------------------------
      SUBROUTINE SUMXPE_2D(arr, arr_master, increment)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:) :: arr
      REAL*8, DIMENSION(:,:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_MPI
         arr_size = size(arr)
         if(increment_) then
            if(am_i_root(COMMUNICATOR)) then
               allocate(arr_tmp(arr_size))
            else
               allocate(arr_tmp(1))
            end if
            call MPI_Reduce(arr,arr_tmp,arr_size, &
    &           MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &           COMMUNICATOR, ierr)
            if(am_i_root(COMMUNICATOR)) then
              arr_master = arr_master + reshape(arr_tmp,shape(arr))
            endif
            deallocate(arr_tmp)
         else
            call MPI_Reduce(arr,arr_master,arr_size, &
    &           MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &           COMMUNICATOR, ierr)
         endif
#else 
         if(increment_) then
            arr_master = arr_master + arr
         else
            arr_master = arr
         endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr is overwritten by itself after reduction
#ifdef USE_MPI
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call MPI_Reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        COMMUNICATOR, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_2D

! ----------------------------------------------------------------------
      SUBROUTINE SUMXPE_3D(arr, arr_master, increment)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:,:) :: arr
      REAL*8, DIMENSION(:,:,:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_MPI
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root(COMMUNICATOR)) then
           allocate(arr_tmp(arr_size))
        else
           allocate(arr_tmp(1))
        end if
        call MPI_Reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,COMMUNICATOR, ierr)
        if(am_i_root(COMMUNICATOR)) then
          arr_master = arr_master + reshape(arr_tmp,shape(arr))
        endif
        deallocate(arr_tmp)
      else
        call MPI_Reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,COMMUNICATOR, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr  is overwritten by itself after reduction
#ifdef USE_MPI
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call MPI_Reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        COMMUNICATOR, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_3D

! ----------------------------------------------------------------------
      SUBROUTINE SUMXPE_4D(arr, arr_master, increment)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:,:,:) :: arr
      REAL*8, DIMENSION(:,:,:,:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_MPI
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root(COMMUNICATOR)) then
           allocate(arr_tmp(arr_size))
        else
           allocate(arr_tmp(1))
        end if
        call MPI_Reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,COMMUNICATOR, ierr)
        if(am_i_root(COMMUNICATOR)) then
          arr_master = arr_master + reshape(arr_tmp,shape(arr))
        endif
        deallocate(arr_tmp)
      else
        call MPI_Reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,COMMUNICATOR, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr  is overwritten by itself after reduction
#ifdef USE_MPI
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call MPI_Reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        COMMUNICATOR, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_4D

#ifdef USE_MPI

      subroutine gridPELayout  (GRID, NX, NY)
        type (Dist_Grid), intent(IN) :: grid
        integer, intent(OUT)         :: NX, NY

        NX = 1
        NY = getNumProcesses(grid)

      end subroutine gridPELayout

      subroutine gridRootPELocation  (GRID,NX0,NY0)
        type (Dist_Grid), intent(IN) :: grid
        integer, intent(OUT)          :: NX0, NY0

        NX0 = 0
        NY0 = grid%rank

      end subroutine gridRootPELocation

#endif

! ----------------------------------------------------------------------
      SUBROUTINE HERE(file, line)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      Character(Len=*) :: file
      Integer :: line

      INTEGER :: ierr
      Integer, Allocatable :: lines(:)

#ifdef DEBUG_DECOMP
      CALL LOG_PARALLEL(grid, file, line)
      If (AM_I_ROOT(COMMUNICATOR)) Then
         WRITE(CHECKSUM_UNIT,*)'HERE: ',file, line
         CALL SYS_FLUSH(CHECKSUM_UNIT)
       End If
#ifdef USE_MPI
       ALLOCATE(lines(npes_world))
       Call MPI_Allgather(line, 1, MPI_INTEGER, lines, 1, MPI_INTEGER, &
    &      COMMUNICATOR, ierr)
       If (Any(lines /= line)) &
    &      call stop_model('HERE: synchronization error -severe.',255)
       Deallocate(lines)
#endif
#endif

      CALL barrier()

      END SUBROUTINE HERE

! ----------------------------------------------------------------------
      SUBROUTINE LOG_PARALLEL(distGrid, file, line, i0, i1, x0, x1)
! ----------------------------------------------------------------------
      Use FILEMANAGER, only : nameunit
      IMPLICIT NONE

      TYPE(DIST_GRID), INTENT(IN) :: distGrid
      CHARACTER(Len=*) :: file
      INTEGER          :: line

      INTEGER, OPTIONAL :: i0, i1(:)
      REAL*8, OPTIONAL  :: x0, x1(:)

      INTEGER :: iu
      INTEGER :: n

#ifdef DEBUG_DECOMP
      iu = distGrid%private%log_unit
      WRITE(iu, *) file, line
      If (PRESENT(i0)) WRITE(iu, *) '   i0=',i0
      If (PRESENT(x0)) WRITE(iu, *) '   x0=',x0
      IF (PRESENT(i1)) THEN
         DO n = 1, Size(i1)
            WRITE(iu, '(10x,i4,1x,a,i6)')n, '   i1=',i1(n)
         END DO
      END IF
      IF (PRESENT(x1)) THEN
         DO n = 1, Size(x1)
            WRITE(iu, '(10x,i4,1x,a,e22.16)')n, '   x1=',x1(n)
         END DO
      END IF
      CALL SYS_FLUSH(iu)
#endif

      END SUBROUTINE LOG_PARALLEL

! ----------------------------------------------------------------------
      SUBROUTINE GLOBALMIN_R(distGrid, val, val_min)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: distGrid
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_min

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_min, 1, MPI_DOUBLE_PRECISION,MPI_MIN, &
    &     getMpiCommunicator(DISTGRID), ierr)
#else
      val_min = val
#endif

      END SUBROUTINE

! ----------------------------------------------------------------------
      SUBROUTINE GLOBALMAX_R(distGrid, val, val_max)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: distGrid
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_max

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_max, 1, MPI_DOUBLE_PRECISION,MPI_MAX, &
    &     getMpiCommunicator(distGrid), ierr)
#else
      val_max = val
#endif

      END SUBROUTINE

! ----------------------------------------------------------------------
      SUBROUTINE GLOBALMIN_I(distGrid, val, val_min)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: distGrid
      INTEGER,            INTENT(IN)  :: val
      INTEGER,            INTENT(OUT) :: val_min

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_min, 1, MPI_INTEGER, MPI_MIN, &
    &     getMpiCommunicator(distGrid), ierr)
#else
      val_min = val
#endif

      END SUBROUTINE

! ----------------------------------------------------------------------
      SUBROUTINE GLOBALMAX_I(distGrid, val, val_max)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: distGrid
      INTEGER,            INTENT(IN)  :: val
      INTEGER,            INTENT(OUT) :: val_max

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_max, 1, MPI_INTEGER, MPI_MAX, &
    &     getMpiCommunicator(distGrid), ierr)
#else
      val_max = val
#endif

      END SUBROUTINE

! ----------------------------------------------------------------------
      SUBROUTINE GLOBALMAX_I_1D(distGrid, val, val_max)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: distGrid
      INTEGER,            INTENT(IN)  :: val(:)
      INTEGER,            INTENT(OUT) :: val_max(:)

      INTEGER  :: n,ierr

#ifdef USE_MPI
      n = size(val)
      CALL MPI_Allreduce(val, val_max, n, MPI_INTEGER, MPI_MAX, &
    &     getMpiCommunicator(distGrid), ierr)
#else
      val_max(:) = val(:)
#endif

      END SUBROUTINE


! ----------------------------------------------------------------------
      SUBROUTINE INIT_BAND_PACK_TYPE(grd_src, grd_dst, band_j0,band_j1, &
    &     bandpack)
! initialize the bandpack derived type with the information needed
! for the band_pack procedure to fill output arrays with data from
! J indices band_j0 to band_j1
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_src,grd_dst
      INTEGER, INTENT(IN) :: band_j0,band_j1
      TYPE (BAND_PACK_TYPE), intent(OUT) :: bandpack
#ifdef USE_MPI
      integer, dimension(0:grd_src%npes_world-1) :: &
    &     j0_have,j1_have,j0_requested,j1_requested
      integer :: p, ierr, im,jm, j0send,j1send,j0recv,j1recv,npes

!
! Store some general MPI info (NOTE: for now we assume that the process
! set of grd_src is a subset of that in grd_dst, or vice versa, which
! allows us to take this info from the larger set).
!
      if(getNumProcesses(grd_src) > getNumProcesses(grd_dst)) then
        npes = getNumProcesses(grd_src)
        bandpack%mpi_comm = getMpiCommunicator(grd_src)
      else
        npes = getNumProcesses(grd_dst)
        bandpack%mpi_comm = getMpiCommunicator(grd_dst)
      endif
      bandpack%npes_comm = npes
      allocate(bandpack%j0_send(0:npes-1), &
    &         bandpack%j1_send(0:npes-1), &
    &         bandpack%j0_recv(0:npes-1), &
    &         bandpack%j1_recv(0:npes-1), &
    &         bandpack%scnts(0:npes-1), &
    &         bandpack%sdspl(0:npes-1), &
    &         bandpack%rcnts(0:npes-1), &
    &         bandpack%rdspl(0:npes-1), &
    &         bandpack%sdspl_inplace(0:npes-1), &
    &         bandpack%rdspl_inplace(0:npes-1) &
    &     )
#endif
      bandpack%j_strt = grd_src%j_strt
      bandpack%j_stop = grd_src%j_stop
      bandpack%j_strt_halo = grd_src%j_strt_halo
      bandpack%j_stop_halo = grd_src%j_stop_halo
      bandpack%jband_strt = band_j0
      bandpack%jband_stop = band_j1
#ifdef USE_MPI
      im = 1
      jm = grd_src%jm_world
!
! Set up the MPI send/receive information
!
      call mpi_allgather(grd_src%j_strt,1,MPI_INTEGER,j0_have,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      call mpi_allgather(grd_src%j_stop,1,MPI_INTEGER,j1_have,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      call mpi_allgather(band_j0,1,MPI_INTEGER,j0_requested,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      call mpi_allgather(band_j1,1,MPI_INTEGER,j1_requested,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      do p=0,npes-1
        j0send = max(grd_src%j_strt,j0_requested(p))
        j1send = min(grd_src%j_stop,j1_requested(p))
        bandpack%j0_send(p) = j0send
        bandpack%j1_send(p) = j1send
        if(j0send <= j1send) then
          bandpack%scnts(p) = im*(j1send-j0send+1)
          bandpack%sdspl_inplace(p) = im*(j0send-grd_src%j_strt_halo)
        else
          bandpack%scnts(p) = 0
          bandpack%sdspl_inplace(p) = 0
        endif
        j0recv = max(j0_have(p),band_j0)
        j1recv = min(j1_have(p),band_j1)
        bandpack%j0_recv(p) = j0recv
        bandpack%j1_recv(p) = j1recv
        if(j0recv <= j1recv) then
          bandpack%rcnts(p) = im*(j1recv-j0recv+1)
          bandpack%rdspl_inplace(p) = im*(j0recv-band_j0)
        else
          bandpack%rcnts(p) = 0
          bandpack%rdspl_inplace(p) = 0
        endif
      enddo
      bandpack%rdspl(0) = 0
      bandpack%sdspl(0) = 0
      do p=1,npes-1
        bandpack%sdspl(p) = bandpack%sdspl(p-1)+bandpack%scnts(p-1)
        bandpack%rdspl(p) = bandpack%rdspl(p-1)+bandpack%rcnts(p-1)
      enddo
#endif
      RETURN
      END SUBROUTINE INIT_BAND_PACK_TYPE

! ----------------------------------------------------------------------
      SUBROUTINE BAND_PACK_ij(bandpack,ARR,ARR_band)
!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
! ----------------------------------------------------------------------
      IMPLICIT NONE 
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,bandpack%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,bandpack%jband_strt:)
#ifdef USE_MPI
      integer, allocatable, dimension(:) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr 
      integer :: im,npes,npes_world
      call mpi_comm_size(bandpack%mpi_comm,npes_world,ierr)
      allocate(scnts(0:npes_world-1))
      allocate(sdspl(0:npes_world-1))
      allocate(rcnts(0:npes_world-1))
      allocate(rdspl(0:npes_world-1))
      npes = bandpack%npes_comm
      im = size(arr,1)
      scnts(0:npes-1) = im*bandpack%scnts
      sdspl(0:npes-1) = im*bandpack%sdspl_inplace
      rcnts(0:npes-1) = im*bandpack%rcnts
      rdspl(0:npes-1) = im*bandpack%rdspl_inplace
      call mpi_alltoallv(arr, scnts, sdspl, mpi_double_precision, &
    &                   arr_band, rcnts, rdspl, mpi_double_precision, &
    &                   bandpack%mpi_comm, ierr)
      deallocate(scnts,sdspl,rcnts,rdspl)
#else
      arr_band(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP) = &
    &     arr(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_ij

! ----------------------------------------------------------------------
      SUBROUTINE BAND_PACK_ijl(bandpack,ARR,ARR_band)
!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,bandpack%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,bandpack%jband_strt:,:)
#ifdef USE_MPI
      integer, allocatable, dimension(:) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr,im,lm,i,j,l,n,p
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      integer :: npes,npes_world
      call mpi_comm_size(bandpack%mpi_comm,npes_world,ierr)
      allocate(scnts(0:npes_world-1))
      allocate(sdspl(0:npes_world-1))
      allocate(rcnts(0:npes_world-1))
      allocate(rdspl(0:npes_world-1))
      npes = bandpack%npes_comm
      im = size(arr,1)
      lm = size(arr,3)
      scnts = im*lm*bandpack%scnts
      sdspl = im*lm*bandpack%sdspl
      rcnts = im*lm*bandpack%rcnts
      rdspl = im*lm*bandpack%rdspl
      allocate(bufsend(sum(scnts)),bufrecv(sum(rcnts)))
      n = 0
      do p=0,npes-1
        do l=1,lm
          do j=bandpack%j0_send(p),bandpack%j1_send(p)
            do i=1,im
              n = n + 1
              bufsend(n) = arr(i,j,l)
            enddo
          enddo
        enddo
      enddo
      call mpi_alltoallv(bufsend, scnts, sdspl, mpi_double_precision, &
    &                   bufrecv, rcnts, rdspl, mpi_double_precision, &
    &                   bandpack%mpi_comm, ierr)
      n = 0
      do p=0,npes-1
        do l=1,lm
          do j=bandpack%j0_recv(p),bandpack%j1_recv(p)
            do i=1,im
              n = n + 1
              arr_band(i,j,l) = bufrecv(n)
            enddo
          enddo
        enddo
      enddo
      deallocate(bufsend,bufrecv)
      deallocate(scnts,sdspl,rcnts,rdspl)
#else
      arr_band(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP,:) = &
    &     arr(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP,:)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_ijl

! ----------------------------------------------------------------------
      SUBROUTINE BAND_PACK_COLUMN(bandpack,ARR,ARR_band)
!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,:,bandpack%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,:,bandpack%jband_strt:)
#ifdef USE_MPI
      integer, allocatable, dimension(:) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr,im,lm
      integer :: npes,npes_world
      call mpi_comm_size(bandpack%mpi_comm,npes_world,ierr)
      allocate(scnts(0:npes_world-1))
      allocate(sdspl(0:npes_world-1))
      allocate(rcnts(0:npes_world-1))
      allocate(rdspl(0:npes_world-1))
      npes = bandpack%npes_comm
      im = size(arr,2)
      lm = size(arr,1)
      scnts(0:npes-1) = im*lm*bandpack%scnts
      sdspl(0:npes-1) = im*lm*bandpack%sdspl_inplace
      rcnts(0:npes-1) = im*lm*bandpack%rcnts
      rdspl(0:npes-1) = im*lm*bandpack%rdspl_inplace
      call mpi_alltoallv(arr, scnts, sdspl, mpi_double_precision, &
    &                   arr_band, rcnts, rdspl, mpi_double_precision, &
    &                   bandpack%mpi_comm, ierr)
      deallocate(scnts,sdspl,rcnts,rdspl)
#else
      arr_band(:,:,bandpack%JBAND_STRT:bandpack%JBAND_STOP) = &
    &     arr(:,:,bandpack%JBAND_STRT:bandpack%JBAND_STOP)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_COLUMN

! ----------------------------------------------------------------------
      SUBROUTINE broadcast_0D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(InOut) :: arr

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,1,MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif

      END SUBROUTINE broadcast_0D

! ----------------------------------------------------------------------
      SUBROUTINE broadcast_1D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(InOut) :: arr(:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif

      END SUBROUTINE broadcast_1D

! ----------------------------------------------------------------------
      SUBROUTINE broadcast_2D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(InOut) :: arr(:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif

      END SUBROUTINE broadcast_2D

! ----------------------------------------------------------------------
      SUBROUTINE broadcast_3D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(InOut) :: arr(:,:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif

      END SUBROUTINE broadcast_3D

! ----------------------------------------------------------------------
      SUBROUTINE broadcast_4D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(InOut) :: arr(:,:,:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif

      END SUBROUTINE broadcast_4D

! ----------------------------------------------------------------------
      SUBROUTINE ibroadcast_0D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(InOut) :: arr
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,1,MPI_INTEGER,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif
      END SUBROUTINE ibroadcast_0D

! ----------------------------------------------------------------------
      SUBROUTINE ibroadcast_0D_world(arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      Integer, Intent(InOut) :: arr
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,1,MPI_INTEGER,root, &
    &     COMMUNICATOR, ierr)
#endif
      END SUBROUTINE ibroadcast_0D_world

! ----------------------------------------------------------------------
      SUBROUTINE ibroadcast_1D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(InOut) :: arr(:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER, root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif
      END SUBROUTINE ibroadcast_1D

! ----------------------------------------------------------------------
      SUBROUTINE ibroadcast_2D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(InOut) :: arr(:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER ,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif
      END SUBROUTINE ibroadcast_2D

! ----------------------------------------------------------------------
      SUBROUTINE ibroadcast_3D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(InOut) :: arr(:,:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER, root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif
      END SUBROUTINE ibroadcast_3D

! ----------------------------------------------------------------------
      SUBROUTINE ibroadcast_4D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(InOut) :: arr(:,:,:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER ,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif
      END SUBROUTINE ibroadcast_4D


! ----------------------------------------------------------------------
      SUBROUTINE TRANSPOSE_ijk(grid, x_in, x_out, reverse)
! ----------------------------------------------------------------------
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x_in(:,grid%J_STRT_HALO:,:)
      REAL*8 :: x_out(:,:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER, DIMENSION(0:grid%NPES_WORLD-1) :: I0, J0, I1, J1
      INTEGER, DIMENSION(0:GRID%NPES_WORLD-1) :: scnts, rcnts, sdspl, rdspl
      REAL*8, ALLOCATABLE :: sbuf(:), rbuf(:)
#ifdef USE_MPI
      TYPE (AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J,nk,k
      INTEGER :: ierr, pc, p
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt, npes
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_MPI
      If (reverse_) Then
         X_IN(:,1:grid%JM_WORLD,:) = X_OUT(:,1:grid%JM_WORLD,:)
      Else
         X_OUT(:,1:grid%JM_WORLD,:) = X_IN(:,1:grid%JM_WORLD,:)
      End If
#else
      npes = getNumProcesses(grid)

      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:grid%npes_world-1,3))
      call getAxisIndex(grid, AI)

      DO p = 0, npes - 1
         J0(p) = AI(p,2)%min
         J1(p) = AI(p,2)%max
      END DO
      DEALLOCATE(AI)

      ni_loc = I1(grid%rank) - I0(grid%rank) + 1
      nj_loc = J1(grid%rank) - J0(grid%rank) + 1

      nk = SIZE(X_IN,3)

      ALLOCATE(rbuf(grid%JM_WORLD * ni_loc * nk))
      ALLOCATE(sbuf(grid%IM_WORLD * nj_loc * nk))

      sdspl(0) = 0
      rdspl(0) = 0
      icnt = 0
      DO p = 0, npes -1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  rbuf(icnt) = X_out(i,j,k)
               END DO
            END DO
         ELSE
            DO j = J0(grid%rank), J1(grid%rank)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  sbuf(icnt) = X_in(i,j,k)
               End Do
            END DO
         END IF
         nip = I1(p) - I0(p) + 1
         njp = J1(p) - J0(p) + 1
         scnts(p) = nj_loc * nip * nk
         rcnts(p) = ni_loc * njp * nk
         If (p > 0) sdspl(p) = sdspl(p-1) + scnts(p-1)
         If (p > 0) rdspl(p) = rdspl(p-1) + rcnts(p-1)
       END DO
      END DO

      If (reverse_) Then
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      End If

      icnt = 0
      DO p = 0, npes - 1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(grid%rank), J1(grid%rank)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  X_in(i,j,k) = sbuf(icnt)
               End Do
            END DO
         Else
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  X_out(i,j,k) = rbuf(icnt)
               END DO
            END DO
         End If
       END DO
      END DO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)
#endif

      END SUBROUTINE TRANSPOSE_ijk

! ----------------------------------------------------------------------
      SUBROUTINE TRANSPOSE_ij(grid, x_in, x_out, reverse)
! ----------------------------------------------------------------------
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x_in(:,grid%J_STRT_HALO:)
      REAL*8 :: x_out(:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER, DIMENSION(0:grid%NPES_WORLD-1) :: I0, J0, I1, J1
      INTEGER, DIMENSION(0:GRID%NPES_WORLD-1) :: scnts, rcnts, sdspl, rdspl
      REAL*8, ALLOCATABLE :: sbuf(:), rbuf(:)
#ifdef USE_MPI
      TYPE (AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J
      INTEGER :: ierr, p
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt, npes
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_MPI
      If (reverse_) Then
         X_IN(:,1:grid%JM_WORLD) = X_OUT(:,1:grid%JM_WORLD)
      Else
         X_OUT(:,1:grid%JM_WORLD) = X_IN(:,1:grid%JM_WORLD)
      End If
#else
      npes = getNumProcesses(grid)

      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:grid%npes_world-1,3))
      call getAxisIndex(grid, AI)

      DO p = 0, npes - 1
         J0(p) = AI(p,2)%min
         J1(p) = AI(p,2)%max
      END DO
      DEALLOCATE(AI)

      ni_loc = I1(grid%rank) - I0(grid%rank) + 1
      nj_loc = J1(grid%rank) - J0(grid%rank) + 1


      ALLOCATE(rbuf(grid%JM_WORLD * ni_loc ))
      ALLOCATE(sbuf(grid%IM_WORLD * nj_loc ))

      sdspl(0) = 0
      rdspl(0) = 0
      icnt = 0
      DO p = 0, npes -1
         If (reverse_) Then
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  rbuf(icnt) = X_out(i,j)
               END DO
            END DO
         ELSE
            DO j = J0(grid%rank), J1(grid%rank)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  sbuf(icnt) = X_in(i,j)
               End Do
            END DO
         END IF
         nip = I1(p) - I0(p) + 1
         njp = J1(p) - J0(p) + 1
         scnts(p) = nj_loc * nip
         rcnts(p) = ni_loc * njp
         If (p > 0) sdspl(p) = sdspl(p-1) + scnts(p-1)
         If (p > 0) rdspl(p) = rdspl(p-1) + rcnts(p-1)
      END DO

      If (reverse_) Then
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      End If

      icnt = 0
      DO p = 0, npes - 1
         If (reverse_) Then
            DO j = J0(grid%rank), J1(grid%rank)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  X_in(i,j) = sbuf(icnt)
               End Do
            END DO
         Else
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  X_out(i,j) = rbuf(icnt)
               END DO
            END DO
         End If
      END DO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)
#endif

      END SUBROUTINE TRANSPOSE_ij

! ----------------------------------------------------------------------
      SUBROUTINE TRANSPOSE_COLUMN(grid, x, x_tr, reverse)
! ----------------------------------------------------------------------
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x(:,:,grid%J_STRT_HALO:,:)
      REAL*8 :: x_tr(:,:,:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER, DIMENSION(0:grid%NPES_WORLD-1) :: I0, J0, I1, J1
      INTEGER, DIMENSION(0:GRID%NPES_WORLD-1) :: scnts, rcnts, sdspl, rdspl
      REAL*8, ALLOCATABLE :: sbuf(:,:), rbuf(:,:)
#ifdef USE_MPI
      TYPE (AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J,k
      INTEGER :: ierr, p
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt, npes
      INTEGER :: n, nk
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_MPI
      If (reverse_) Then
         X(:,:,1:grid%JM_WORLD,:) = X_TR(:,:,1:grid%JM_WORLD,:)
      Else
         X_TR(:,:,1:grid%JM_WORLD,:) = X(:,:,1:grid%JM_WORLD,:)
      End If
#else
      npes = getNumProcesses(grid)

      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:grid%npes_world-1,3))
      call getAxisIndex(grid, AI)

      DO p = 0, npes - 1
         J0(p) = AI(p,2)%min
         J1(p) = AI(p,2)%max
      END DO
      DEALLOCATE(AI)

      ni_loc = I1(grid%rank) - I0(grid%rank) + 1
      nj_loc = J1(grid%rank) - J0(grid%rank) + 1

      n  = SIZE(X, 1)
      nk = SIZE(X,4)

      ALLOCATE(rbuf(n, grid%JM_WORLD * ni_loc * nk))
      ALLOCATE(sbuf(n, grid%IM_WORLD * nj_loc * nk))

      sdspl(0) = 0
      rdspl(0) = 0
      icnt = 0
      DO p = 0, npes -1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  rbuf(:,icnt) = X_tr(:,i,j,k)
               END DO
            END DO
         ELSE
            DO j = J0(grid%rank), J1(grid%rank)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  sbuf(:,icnt) = X(:,i,j,k)
               End Do
            END DO
         END IF
         nip = I1(p) - I0(p) + 1
         njp = J1(p) - J0(p) + 1
         scnts(p) = n * nj_loc * nip * nk
         rcnts(p) = n * ni_loc * njp * nk
         If (p > 0) sdspl(p) = sdspl(p-1) + scnts(p-1)
         If (p > 0) rdspl(p) = rdspl(p-1) + rcnts(p-1)
       END DO
      END DO

      If (reverse_) Then
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      End If

      icnt = 0
      DO p = 0, npes - 1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(grid%rank), J1(grid%rank)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  X(:,i,j,k) = sbuf(:,icnt)
               End Do
            END DO
         Else
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  X_tr(:,i,j,k) = rbuf(:,icnt)
               END DO
            END DO
         End If
       END DO
      END DO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)
#endif

      END SUBROUTINE TRANSPOSE_COLUMN

! ----------------------------------------------------------------------
      logical function haveLatitude(distGrid, j)
! ----------------------------------------------------------------------
        type (DIST_GRID), intent(in) :: distGrid
        integer, intent(in) :: j

        haveLatitude = (j >= distGrid%j_strt .and. j <= distGrid%J_STOP)

      end function haveLatitude

! ----------------------------------------------------------------------
      logical function isInLocalSubdomain(distGrid, i, j)
! ----------------------------------------------------------------------
        type (DIST_GRID), intent(in) :: distGrid
        integer, intent(in) :: i
        integer, intent(in) :: j

        isInLocalSubdomain = (i >= distGrid%i_strt .and. i <= distGrid%i_stop)
        isInLocalSubdomain = isInLocalSubdomain .and. &
             (j >= distGrid%j_strt .and. j <= distGrid%j_stop)

      end function isInLocalSubdomain

! ----------------------------------------------------------------------
      subroutine SEND_TO_J_1D(distGrid, arr, j_dest, tag)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(In) :: arr(:)
      Integer, Intent(In) :: j_dest, tag
      INTEGER :: ierr
#ifdef USE_MPI
      call MPI_Send(arr, Size(arr), MPI_DOUBLE_PRECISION, &
    &     distGrid%private%lookup_pet(j_dest), tag, COMMUNICATOR, ierr)
#endif
      end subroutine SEND_TO_J_1D

! ----------------------------------------------------------------------
      subroutine ISEND_TO_J_0D(distGrid, arr, j_dest, tag)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(In) :: arr
      Integer, Intent(In) :: j_dest, tag
      INTEGER :: ierr
#ifdef USE_MPI
      call MPI_Send(arr, 1, MPI_INTEGER, &
    &     distGrid%private%lookup_pet(j_dest), tag, COMMUNICATOR, ierr)
#endif
      end subroutine ISEND_TO_J_0D

! ----------------------------------------------------------------------
      subroutine RECV_FROM_J_1D(distGrid, arr, j_src, tag)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(In) :: arr(:)
      Integer, Intent(In) :: j_src, tag
#ifdef USE_MPI
      INTEGER :: ierr, status(MPI_STATUS_SIZE)
      call MPI_Recv(arr, Size(arr), MPI_DOUBLE_PRECISION, &
    &     distGrid%private%lookup_pet(j_src), tag, COMMUNICATOR, status, ierr)
#endif
      end subroutine RECV_FROM_J_1D

! ----------------------------------------------------------------------
      subroutine IRECV_FROM_J_0D(distGrid, arr, j_src, tag)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(In) :: arr
      Integer, Intent(In) :: j_src, tag
#ifdef USE_MPI
      INTEGER :: ierr, status(MPI_STATUS_SIZE)
      call MPI_Recv(arr, 1, MPI_INTEGER, &
    &     distGrid%private%lookup_pet(j_src), tag, COMMUNICATOR, status, ierr)
#endif
      end subroutine IRECV_FROM_J_0D

! ----------------------------------------------------------------------
! Helper function to handle optional arguments related to periodic boundaries
      logical function isPeriodic(override)
! ----------------------------------------------------------------------
        logical, optional, intent(in) :: override

        isPeriodic = .false.
        if (present(override)) isPeriodic = override

      end function isPeriodic

! ----------------------------------------------------------------------
      subroutine barrier()
! ----------------------------------------------------------------------
        integer :: rc
#ifdef USE_ESMF
        call ESMF_VMbarrier(modelE_vm, rc=rc)
#else
#ifdef USE_MPI
        call MPI_Barrier(communicator, rc)
#endif
#endif
      end subroutine barrier

! ----------------------------------------------------------------------
      subroutine setMpiCommunicator(this, comm)
! ----------------------------------------------------------------------
        type (dist_grid), intent(inout) :: this
        integer, intent(in) :: comm
        this%private%MPI_COMM = comm
      end subroutine setMpiCommunicator

! ----------------------------------------------------------------------
      subroutine setCommunicator(comm)
! ----------------------------------------------------------------------
        integer, intent(in) :: comm
        communicator = comm
      end subroutine setCommunicator

! ----------------------------------------------------------------------
      integer function getMpiCommunicator(this) result(comm)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        comm = this%private%mpi_comm
      end function getMpiCommunicator

! ----------------------------------------------------------------------
      integer function getNumProcesses(this) result(numProcesses)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        numProcesses = this%private%numProcesses
      end function getNumProcesses

! ----------------------------------------------------------------------
      integer function getNumAllProcesses(this) result(numAllProcesses)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        numAllProcesses = this%private%numAllProcesses
      end function getNumAllProcesses


! ----------------------------------------------------------------------
      subroutine incrementMpiTag(this)
! ----------------------------------------------------------------------
        type (dist_grid), intent(inout) :: this
        integer, parameter :: MIN_TAG = 10
        integer, parameter :: MAX_TAG = 128

        integer :: tag
        tag = this%private%mpi_tag
        this%private%mpi_tag = max(mod(tag,MAX_TAG),MIN_TAG) + 1
      end subroutine incrementMpiTag

! ----------------------------------------------------------------------
      integer function getMpiTag(this) result(mpiTag)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        mpiTag = this%private%mpi_tag
      end function getMpiTag

! ----------------------------------------------------------------------
      logical function hasSouthPole(this)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        hasSouthPole = this%private%hasSouthPole
      end function hasSouthPole

! ----------------------------------------------------------------------
      logical function hasNorthPole(this)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        hasNorthPole = this%private%hasNorthPole
      end function hasNorthPole

! ----------------------------------------------------------------------
      logical function hasPeriodicBC(this)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        hasPeriodicBC = this%private%periodicBC
      end function hasPeriodicBC

! ----------------------------------------------------------------------
      integer function getLogUnit()
! ----------------------------------------------------------------------
        use FileManager, only: openUnit
        character(len=40) :: logFileName

        integer, parameter :: UNINITIALIZED = -1
        integer, save :: logUnit = UNINITIALIZED

        if (logUnit == UNINITIALIZED) then
          write(logFileName,'(a,i4.4)') 'debug.', rank
          call openUnit(logFileName, logUnit, qbin=.false., qold=.false.)
        end if

        getLogUnit = logUnit

      end function getLogUnit

!-------------------------------------------------------------------------------
  Subroutine abort(line,rc)
!-------------------------------------------------------------------------------
    Implicit None
    Integer, Intent(In) :: line
    Integer, Intent(In) :: rc

    Character(len=100) :: buf

    Write(buf,*)__FILE__,' failure at line',line,'\n error code:', rc
    Call stop_model(Trim(buf),99)

  End Subroutine abort

END MODULE dist_grid_mod

