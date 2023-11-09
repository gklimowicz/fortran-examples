module hybrid_mpi_omp_coupler
  ! Variables/entities used from modelE

  USE SEAICE, only:     fsss   ! scalar
  USE SEAICE, only:     tfrez  ! procedure

  USE MODEL_COM, only: dtsrc,itime,iyear1,nday,jdendofm,aMON,modelEclock

  USE GEOM, only: dxyp
  USE CONSTANT,  only: lhm,shi,shw

  implicit none
  private

  public :: init_hybrid_coupler
  public :: gatherDistributedQuantities
  public :: scatterDistributedQuantities
  public :: gatherDistributedQuantities_IO
  public :: scatterDistributedQuantities_IO

  ! Global (gathered) quantities related to modelE local (scattered quantities)
  
  ! From FLUXES
  PUBLIC :: e0    
  PUBLIC :: prec  
  PUBLIC :: eprec 
  PUBLIC :: evapor
  PUBLIC :: flowo 
  PUBLIC :: eflowo
  PUBLIC :: dmua  
  PUBLIC :: dmva   
                                                    
  PUBLIC :: erunosi
  PUBLIC :: runosi 
  PUBLIC :: srunosi
  PUBLIC :: runpsi 
  PUBLIC :: srunpsi
  PUBLIC :: dmui   
  PUBLIC :: dmvi   
  PUBLIC :: dmsi   
  PUBLIC :: dhsi   
  PUBLIC :: dssi   
                                                    
  PUBLIC :: gtemp
  PUBLIC :: sss  
  PUBLIC :: mlhc 
  PUBLIC :: ogeoza
  PUBLIC :: uosurf
  PUBLIC :: vosurf

  PUBLIC :: MELTI
  PUBLIC :: EMELTI
  PUBLIC :: SMELTI

  PUBLIC :: gmelt
  PUBLIC :: egmelt
  PUBLIC :: solar

  ! from SEAICE_COM
  PUBLIC :: rsi
  PUBLIC :: msi
  ! from SEAICE 
  PUBLIC :: fsss
  PUBLIC :: tfrez

  ! from MODEL_COM
  PUBLIC :: focean
  PUBLIC :: dtsrc,itime,iyear1,nday,jdendofm,jyear,jmon,jday,jdate,jhour,aMON

  ! from GEOM
  PUBLIC :: dxyp
  ! from CONSTANT
  PUBLIC :: lhm,shi,shw

  real*8, allocatable :: e0(:,:,:)
  real*8, allocatable :: prec(:,:)
  real*8, allocatable :: eprec(:,:)
  real*8, allocatable :: evapor(:,:,:)
  real*8, allocatable :: flowo(:,:)
  real*8, allocatable :: eflowo(:,:)
  real*8, allocatable :: dmua(:,:,:)
  real*8, allocatable :: dmva(:,:,:)

  real*8, allocatable :: erunosi(:,:)
  real*8, allocatable :: runosi(:,:)
  real*8, allocatable :: srunosi(:,:)
  real*8, allocatable :: runpsi(:,:)
  real*8, allocatable :: srunpsi(:,:)
  real*8, allocatable :: dmui(:,:)
  real*8, allocatable :: dmvi(:,:)
  real*8, allocatable :: dmsi(:,:,:)
  real*8, allocatable :: dhsi(:,:,:)
  real*8, allocatable :: dssi(:,:,:)

  real*8, allocatable :: gtemp(:,:,:,:)
  real*8, allocatable :: sss(:,:)
  real*8, allocatable :: mlhc(:,:)
  real*8, allocatable :: ogeoza(:,:)
  real*8, allocatable :: uosurf(:,:)
  real*8, allocatable :: vosurf(:,:)

  real*8, allocatable :: melti(:,:)
  real*8, allocatable :: emelti(:,:)
  real*8, allocatable :: smelti(:,:)

  real*8, allocatable :: gmelt(:,:)
  real*8, allocatable :: egmelt(:,:)
  real*8, allocatable :: solar(:,:,:)

  real*8, allocatable :: rsi(:,:)
  real*8, allocatable :: msi(:,:)

  real*8, allocatable :: focean(:,:)

  ! private
  integer, parameter :: SINGLE_THREADED = 1
  integer :: NUM_THREADS_HYCOM

contains

  subroutine gatherDistributedQuantities()
    use MODEL_COM, only: IM, JM
    USE DOMAIN_DECOMP_1D, ONLY: GRID, PACK_DATA, PACK_COLUMN, PACK_BLOCK
    use hybrid_mpi_omp_renamer
    ! initial 1-processor implementation

    call pack_data(grid,  e0_loc,       e0)
    call pack_data(grid,  prec_loc,     prec)
    call pack_data(grid,  eprec_loc,    eprec)
    call pack_data(grid,  evapor_loc,   evapor)
    call pack_data(grid,  flowo_loc ,   flowo)
    call pack_data(grid,  eflowo_loc,   eflowo)
    call pack_data(grid,  dmua_loc  ,   dmua)
    call pack_data(grid,  dmva_loc  ,   dmva)
                                                                     
    call pack_data(grid, erunosi_loc, erunosi) 
    call pack_data(grid, runosi_loc , runosi )
    call pack_data(grid, srunosi_loc, srunosi) 
    call pack_data(grid, runpsi_loc , runpsi )
    call pack_data(grid, srunpsi_loc, srunpsi)
    call pack_data(grid, dmui_loc   , dmui   )
    call pack_data(grid, dmvi_loc   , dmvi   ) 
    call pack_column(grid, dmsi_loc   , dmsi   ) 
    call pack_column(grid, dhsi_loc   , dhsi   )
    call pack_column(grid, dssi_loc   , dssi   )
                                                                    
    call pack_block(grid, gtemp_loc  , gtemp )
    call pack_data(grid, sss_loc    , sss    )
    call pack_data(grid, mlhc_loc   , mlhc   )
    call pack_data(grid, ogeoza_loc , ogeoza )
    call pack_data(grid, uosurf_loc , uosurf )
    call pack_data(grid, vosurf_loc , vosurf )

    call pack_data(grid, melti_loc  , melti  )
    call pack_data(grid, emelti_loc , emelti )
    call pack_data(grid, smelti_loc , smelti )
                                                                      
    call pack_data(grid, gmelt_loc,   gmelt)
    call pack_data(grid, egmelt_loc,  egmelt)  
    call pack_column(grid, solar_loc,   solar)   
                                                                     
    call pack_data(grid, rsi_loc,     rsi )   
    call pack_data(grid, msi_loc,     msi )
    
    call pack_data(grid, focean_loc, focean)

  end subroutine gatherDistributedQuantities

  subroutine scatterDistributedQuantities()
    use MODEL_COM, only: IM, JM
    USE DOMAIN_DECOMP_1D, ONLY: GRID, UNPACK_DATA, UNPACK_COLUMN, UNPACK_BLOCK
    use hybrid_mpi_omp_renamer

    ! The following quantities are import only, and do not need to be
    ! scattered back to modelE:
    ! e0, nprec, eprec, evapor, flowo, eflowo, dmua, dmva
    ! erunosi, runosi, runpsi, srunpsi, dmui, dmvi
    ! mlhc
    ! melti, emelti, smelti, gmelt, egmelt, solar
    ! rsi, msi, focean

    call unpack_column(grid, dmsi   ,  dmsi_loc) ! export
    call unpack_column(grid, dhsi   ,  dhsi_loc) ! export
    call unpack_column(grid, dssi   ,  dssi_loc) ! export
                                                                    
    call unpack_block(grid, gtemp ,  gtemp_loc)   ! export
    call unpack_data(grid, sss    ,  sss_loc)     ! export
    call unpack_data(grid, ogeoza ,  ogeoza_loc)
    call unpack_data(grid, uosurf ,  uosurf_loc)  ! export
    call unpack_data(grid, vosurf ,  vosurf_loc)  ! export

  end subroutine scatterDistributedQuantities

  subroutine gatherDistributedQuantities_IO()
    use MODEL_COM, only: IM, JM
    USE DOMAIN_DECOMP_1D, ONLY: GRID, PACK_DATA, PACK_COLUMN, PACK_BLOCK
    use hybrid_mpi_omp_renamer
    ! initial 1-processor implementation

    call pack_column(grid, dmsi_loc   , dmsi   ) 
    call pack_column(grid, dhsi_loc   , dhsi   )
    call pack_column(grid, dssi_loc   , dssi   )
                                                                 
    call pack_block(grid, gtemp_loc  , gtemp )
    call pack_data(grid, sss_loc    , sss    )
    call pack_data(grid, ogeoza_loc , ogeoza )
    call pack_data(grid, uosurf_loc , uosurf )
    call pack_data(grid, vosurf_loc , vosurf )

  end subroutine gatherDistributedQuantities_IO

  subroutine scatterDistributedQuantities_IO()
    use MODEL_COM, only: IM, JM
    USE DOMAIN_DECOMP_1D, ONLY: GRID, UNPACK_DATA, UNPACK_COLUMN
    use hybrid_mpi_omp_renamer
    ! initial 1-processor implementation

    ! list identical with that used to couple run step
    call scatterDistributedQuantities()

  end subroutine scatterDistributedQuantities_IO

  subroutine init_hybrid_coupler()
    use MODEL_COM, only: IM, JM
    use FLUXES, only: NSTYPE
    use domain_decomp_1d, only: AM_I_ROOT

    integer :: OMP_NUM_THREADS
    integer, parameter :: BUF_LEN=20
    integer :: ier
    character(len=BUF_LEN) :: buf

    if (.not. AM_I_ROOT()) then
       return
    end if

    ! must be root process
!$  call PXFGETENV('NUM_THREADS_HYCOM', len('NUM_THREADS_HYCOM'), buf, len(buf), ier)
!$  read(buf,*) NUM_THREADS_HYCOM
!$  print*,'HYCOM using ',NUM_THREADS_HYCOM, ' threads.'
!$  call omp_set_num_threads(NUM_THREADS_HYCOM)

    allocate(e0(IM,JM,NSTYPE))
    allocate(prec(IM,JM))
    allocate(eprec(IM,JM))
    allocate(evapor(IM,JM,NSTYPE))
    allocate(flowo(IM,JM))
    allocate(eflowo(IM,JM))
    allocate(dmua(IM,JM,NSTYPE))
    allocate(dmva(IM,JM,NSTYPE))
  
    allocate(erunosi(IM,JM))
    allocate(runosi(IM,JM))
    allocate(srunosi(IM,JM))
    allocate(runpsi(IM,JM))
    allocate(srunpsi(IM,JM))
    allocate(dmui(IM,JM))
    allocate(dmvi(IM,JM))
    allocate(dmsi(2,IM,JM))
    allocate(dhsi(2,IM,JM))
    allocate(dssi(2,IM,JM))

    allocate(gtemp(2, NSTYPE, IM, JM))
    gtemp=0
    allocate(sss(IM,JM))
    allocate(mlhc(IM,JM))
    allocate(ogeoza(IM,JM))
    allocate(uosurf(IM,JM))
    allocate(vosurf(IM,JM))

    allocate(melti(IM,JM))
    allocate(emelti(IM,JM))
    allocate(smelti(IM,JM))

    allocate(gmelt(IM,JM))
    allocate(egmelt(IM,JM))
    allocate(solar(3,IM,JM))

    allocate(rsi(IM,JM))
    allocate(msi(IM,JM))

    allocate(focean(IM,JM))

  end subroutine init_hybrid_coupler

  subroutine finalize()
    use domain_decomp_1d, only: AM_I_ROOT

    if (.not. AM_I_ROOT()) return

    deallocate(e0)
    deallocate(prec)
    deallocate(eprec)
    deallocate(evapor)
    deallocate(flowo)
    deallocate(eflowo)
    deallocate(dmua)
    deallocate(dmva)
  
    deallocate(erunosi)
    deallocate(runosi)
    deallocate(srunosi)
    deallocate(runpsi)
    deallocate(srunpsi)
    deallocate(dmui)
    deallocate(dmvi)
    deallocate(dmsi)
    deallocate(dhsi)
    deallocate(dssi)

    deallocate(gtemp)
    deallocate(sss)
    deallocate(mlhc)
    deallocate(ogeoza)
    deallocate(uosurf)
    deallocate(vosurf)

    deallocate(melti)
    deallocate(emelti)
    deallocate(smelti)

    deallocate(gmelt)
    deallocate(egmelt)
    deallocate(solar)

    deallocate(rsi)
    deallocate(msi)

    deallocate(focean)

  end subroutine finalize

end module hybrid_mpi_omp_coupler
