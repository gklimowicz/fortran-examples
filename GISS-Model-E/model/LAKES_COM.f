#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif
      MODULE LAKES_COM
!@sum  LAKES_COM model variables for Lake/Rivers module
!@auth Gavin Schmidt
#ifdef TRACERS_WATER
      USE TRACER_COM, only : NTM
#endif
      USE SEAICE_COM, only : icestate,iceocn_xchng_vars
      IMPLICIT NONE
      SAVE
!@var MWL mass of lake water (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: MWL
!@var GML total enthalpy of lake (J)
      REAL*8, ALLOCATABLE,  DIMENSION(:,:) :: GML
!@var TLAKE temperature of lake (C)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: TLAKE
!@var MLDLK mixed layer depth in lake (m)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: MLDLK
!@var FLAKE variable lake fraction (1)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: FLAKE
!@var TANLK tan(alpha) = slope for conical lake (1)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: TANLK
!@var SVFLAKE previous lake fraction (1)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: SVFLAKE

!@var DLAKE depth of lake (m)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: DLAKE
!@var GLAKE like GML but per unit area of lake (J/m2)
      REAL*8,  ALLOCATABLE, DIMENSION(:,:) :: GLAKE
!@var HLAKE lake sill depth (m)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: HLAKE

#ifdef TRACERS_WATER
!@var TRLAKE tracer amount in each lake level (kg)      
Crgr      REAL*8,  ALLOCATABLE, DIMENSION(NTM,2,:,:) :: TRLAKE
      REAL*8,  ALLOCATABLE, DIMENSION(:,:,:,:) :: TRLAKE
#endif

!      type(icestate) :: lki_state ! not yet: lake-private ice variables

!@var icelak derived-type strucure containing variables
!@+   (or pointers thereto) needed for lake-ice interactions.
      type(iceocn_xchng_vars) :: icelak

      target :: flake,mldlk,dlake,glake,tlake

!@param NRVRMX Max No. of named rivers
      INTEGER, PARAMETER :: NRVRMX = 42
!@var NRVR actual No. of named rivers
      INTEGER :: NRVR
!@var IRVRMTH,JRVRMTH indexes for named river mouths
      INTEGER, DIMENSION(NRVRMX) :: IRVRMTH,JRVRMTH
!@var NAMERVR Named rivers
      CHARACTER*8, DIMENSION(NRVRMX) :: NAMERVR
!@var RVROUT Discharges from named rivers
      REAL*8, DIMENSION(NRVRMX) :: RVROUT

      END MODULE LAKES_COM


       SUBROUTINE ALLOC_LAKES_COM (GRID)
C23456789012345678901234567890123456789012345678901234567890123456789012
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth Raul Garza-Robles
      USE DOMAIN_DECOMP_ATM, only: DIST_GRID, getDomainBounds
      USE LAKES_COM, ONLY: MWL, GML, TLAKE, MLDLK, FLAKE, TANLK, SVFLAKE
     &     ,HLAKE,DLAKE,GLAKE,icelak
#ifdef TRACERS_WATER
      USE TRACER_COM, only : NTM
      USE LAKES_COM, ONLY:  TRLAKE
#endif
      USE EXCHANGE_TYPES, only : alloc_xchng_vars
      USE SEAICE_COM, only : alloc_icestate_type,si_atm
      USE FLUXES, only : atmice
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid
      INTEGER :: I_0H,I_1H, J_0H,J_1H
      INTEGER IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE ( MWL(I_0H:I_1H,J_0H:J_1H),
     *           GML(I_0H:I_1H,J_0H:J_1H),
     *           TLAKE(I_0H:I_1H,J_0H:J_1H),
     *           MLDLK(I_0H:I_1H,J_0H:J_1H),
     *           FLAKE(I_0H:I_1H,J_0H:J_1H),
     *           TANLK(I_0H:I_1H,J_0H:J_1H),
     *           SVFLAKE(I_0H:I_1H,J_0H:J_1H),
     &           HLAKE(I_0H:I_1H,J_0H:J_1H),
     &           DLAKE(I_0H:I_1H,J_0H:J_1H),
     &           GLAKE(I_0H:I_1H,J_0H:J_1H),
     *           STAT=IER
     *            )


#ifdef TRACERS_WATER
!@var TRLAKE tracer amount in each lake level (kg)
      ALLOCATE( TRLAKE(NTM,2,I_0H:I_1H,J_0H:J_1H)
     * , STAT=IER)
      si_atm % ntm = ntm
      icelak % ntm = ntm
#endif

      call alloc_icestate_type(grid,si_atm,'LAKES')
      !atmice%rsi   => si_atm%rsi
      !atmice%snowi => si_atm%snowi

      call alloc_xchng_vars(grid,icelak)

      deallocate(icelak%fwater); icelak%fwater  => flake

      icelak%mldlk   => mldlk
      icelak%dlake   => dlake
      icelak%glake   => glake

      RETURN
      END SUBROUTINE ALLOC_LAKES_COM

      subroutine def_rsf_lakes(fid)
!@sum  def_rsf_lakes defines lake array structure in restart files
!@auth M. Kelley
!@ver  beta
      use lakes_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use conserv_diags
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,mldlk,'mldlk(dist_im,dist_jm)')
      call defvar(grid,fid,mwl,'mwl(dist_im,dist_jm)')
      call defvar(grid,fid,tlake,'tlake(dist_im,dist_jm)')
      call defvar(grid,fid,gml,'gml(dist_im,dist_jm)')
      call defvar(grid,fid,flake,'flake(dist_im,dist_jm)')
#ifdef TRACERS_WATER
      call defvar(grid,fid,trlake,'trlake(ntm,d2,dist_im,dist_jm)')
#endif
      call declare_conserv_diags( grid, fid, 'wliql(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'eliql(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'wlaki(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'elaki(dist_im,dist_jm)' )
      return
      end subroutine def_rsf_lakes

      subroutine new_io_lakes(fid,iaction)
!@sum  new_io_lakes read/write lake arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use lakes_com
      use conserv_diags
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      external conserv_LKM, conserv_LKE
      external conserv_LMSI, conserv_LHSI
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 'mldlk', mldlk)
        call write_dist_data(grid, fid, 'mwl',   mwl)
        call write_dist_data(grid, fid, 'tlake', tlake)
        call write_dist_data(grid, fid, 'gml',   gml)
        call write_dist_data(grid, fid, 'flake', flake)
#ifdef TRACERS_WATER
        call write_dist_data(grid, fid, 'trlake', trlake, jdim=4)
#endif
        call dump_conserv_diags( grid, fid, 'wliql', conserv_LKM )
        call dump_conserv_diags( grid, fid, 'eliql', conserv_LKE )
        call dump_conserv_diags( grid, fid, 'wlaki', conserv_LMSI )
        call dump_conserv_diags( grid, fid, 'elaki', conserv_LHSI )
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'mldlk', mldlk)
        call read_dist_data(grid, fid, 'mwl',   mwl)
        call read_dist_data(grid, fid, 'tlake', tlake)
        call read_dist_data(grid, fid, 'gml',   gml)
        call read_dist_data(grid, fid, 'flake', flake)
#ifdef TRACERS_WATER
        call read_dist_data(grid, fid, 'trlake', trlake, jdim=4)
#endif
      end select
      return
      end subroutine new_io_lakes

      subroutine read_agrice_ic
!@sum   read_agrice_ic read atm-grid floating ice initial conditions file.
      use model_com, only : ioread
      use pario, only : par_open,par_close
      use domain_decomp_atm, only : grid
      use filemanager, only : file_exists
      implicit none
      integer :: fid
      if(file_exists('GIC')) then
        fid = par_open(grid,'GIC','read')
      elseif(file_exists('AICEIC')) then
        fid = par_open(grid,'AICEIC','read')
      else
        return
      endif
      call new_io_agrice (fid,ioread)
      call par_close(grid,fid)
      return
      end subroutine read_agrice_ic

      subroutine def_rsf_agrice(fid)
!@sum  def_rsf_lakes defines atm-grid floating ice arrays in restart files
!@auth M. Kelley
!@ver  beta
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use seaice_com, only : si_atm
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,si_atm%rsi,'rsi_atm(dist_im,dist_jm)')
      call defvar(grid,fid,si_atm%snowi,'snowi_atm(dist_im,dist_jm)')
      call defvar(grid,fid,si_atm%msi,'msi_atm(dist_im,dist_jm)')
      call defvar(grid,fid,si_atm%pond_melt,
     &     'pond_melt_atm(dist_im,dist_jm)')
      call defvar(grid,fid,si_atm%flag_dsws,
     &     'flag_dsws_atm(dist_im,dist_jm)')
      call defvar(grid,fid,si_atm%hsi,'hsi_atm(lmi,dist_im,dist_jm)')
      call defvar(grid,fid,si_atm%ssi,'ssi_atm(lmi,dist_im,dist_jm)')
#ifdef TRACERS_WATER
      call defvar(grid,fid,si_atm%trsi,
     &     'trsi_atm(ntm,lmi,dist_im,dist_jm)')
#endif
      return
      end subroutine def_rsf_agrice

      subroutine new_io_agrice(fid,iaction)
!@sum  new_io_lakes read/write atm-grid floating ice arrays from/to rsf
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use model_com, only : ioread,iowrite
      use seaice_com, only : si_atm
      use fluxes, only : atmice
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call seaice_to_atmgrid(atmice)
        call write_dist_data(grid, fid, 'rsi_atm', si_atm%rsi)
        call write_dist_data(grid, fid, 'snowi_atm', si_atm%snowi)
        call write_dist_data(grid, fid, 'msi_atm', si_atm%msi)
        call write_dist_data(grid, fid, 'pond_melt_atm',
     &       si_atm%pond_melt)
        call write_dist_data(grid, fid, 'flag_dsws_atm',
     &       si_atm%flag_dsws)
        call write_dist_data(grid, fid, 'hsi_atm', si_atm%hsi, jdim=3)
        call write_dist_data(grid, fid, 'ssi_atm', si_atm%ssi, jdim=3)
#ifdef TRACERS_WATER
        call write_dist_data(grid, fid, 'trsi_atm', si_atm%trsi, jdim=4)
#endif

      case (ioread)            ! input from restart file
        si_atm%rsi(:,:) = -1d30
        call read_dist_data(grid, fid, 'rsi_atm', si_atm%rsi)
        if(all(si_atm%rsi(:,:)==-1d30)) then
          call stop_model('new_io_agrice: use a restart file '//
     &         'containing atm-grid ice fields',255)
        endif
        si_atm%RSISAVE(:,:)=si_atm%RSI(:,:)
        call read_dist_data(grid, fid, 'snowi_atm', si_atm%snowi)
        call read_dist_data(grid, fid, 'msi_atm', si_atm%msi)
        call read_dist_data(grid, fid, 'pond_melt_atm',
     &       si_atm%pond_melt)
        call read_dist_data(grid, fid, 'flag_dsws_atm',
     &       si_atm%flag_dsws)
        call read_dist_data(grid, fid, 'hsi_atm', si_atm%hsi, jdim=3)
        call read_dist_data(grid, fid, 'ssi_atm', si_atm%ssi, jdim=3)
#ifdef TRACERS_WATER
        call read_dist_data(grid, fid, 'trsi_atm', si_atm%trsi, jdim=4)
#endif
      end select
      return
      end subroutine new_io_agrice

      subroutine def_meta_rvracc(fid)
!@sum  def_meta_rvracc defines river metadata in acc files
!@auth M. Kelley
!@ver  beta
      use lakes_com, only : nrvr,rvrout,namervr
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_attr
      implicit none
      integer :: fid         !@var fid file id
      if(nrvr.lt.1) return
      call defvar(grid,fid,rvrout(1:nrvr),'rvr(nrvr)')
      call write_attr(grid,fid,'rvr','reduction','sum')
      call defvar(grid,fid,namervr(1:nrvr),'namervr(rvr_strlen,nrvr)')
      return
      end subroutine def_meta_rvracc

      subroutine write_meta_rvracc(fid)
!@sum  write_meta_rvracc write river accumulation metadata to file
!@auth M. Kelley
      use lakes_com, only : nrvr,rvrout,namervr
      use domain_decomp_atm, only : grid
      use pario, only : write_data
      implicit none
      integer fid   !@var fid unit number of read/write
      if(nrvr.lt.1) return
      call write_data(grid,fid,'rvr',rvrout(1:nrvr))
      call write_data(grid,fid,'namervr',namervr(1:nrvr))
      return
      end subroutine write_meta_rvracc
