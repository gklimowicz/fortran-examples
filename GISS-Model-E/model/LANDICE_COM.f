#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif

      MODULE LANDICE_COM
!@sum  LANDICE_COM contains the model arrays for land ice
!@auth Gavin Schmidt
!@ver  2010/10/13
!@cont io_landice
      use cdl_mod, only : cdl_type
      use mdiag_com, only : sname_strlen,units_strlen,lname_strlen
      use iso_c_binding
      IMPLICIT NONE
      SAVE

!@var nhc number of height classes
! Moved to module DOMAIN_DECOMP_ATM
!      type(c_ptr) :: glint2   ! Handle to ice sheet coupler API
      integer :: nhc=1
      REAL*8 :: HC_T_LAPSE_RATE = .008		! Lapse rate to use in T downscaling, K/m

!@usedhp Integer-boolean array that tells whether height points are enabled
!       in each grid cell.  (Generally they're only enabled for grid cells
!       that overlap hi-res ice models)
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: usedhp
!@fhc fraction of landice area in each height class (static for testing purposes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: fhc
#ifdef GLINT2
!@fhc fraction of landice area in each height point (approximate)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: fhc_approx
#endif

!@var ELEVHP: surface elevation, per height class (m)
! ZATMO should be kept consistent with this.
! The value of this ONLY MATTERS for grid cells with landice
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: ELEVHP

!@var SNOWLI snow amount on land ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SNOWLI
!@var TLANDI temperature of each land ice layer (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TLANDI
!@var MDWNIMP downward implicit ice amount accumulator (kg)
!@var EDWNIMP downward implicit energy amount accumulator (J)
!@var FSHGLM = fraction of SH GMELT water; Sum[FSHGLM(:,:)] = 1
!@var FNHGLM = fraction of NH GMELT water; Sum[FNHGLM(:,:)] = 1
      Real*8,Allocatable,Dimension(:,:) ::
     *  MDWNIMP,EDWNIMP, FSHGLM,FNHGLM

#ifdef TRACERS_WATER
!@var TRSNOWLI tracer amount in land ice snow (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRSNOWLI
!@var TRLNDI tracer amount in land ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRLNDI
!@var TDWNIMP downward implicit tracer amount accumulator (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRDWNIMP
#endif

!@param kijhc number of ijhc accumulations
      integer, parameter :: kijhc=16
!@var ijhc accumulations for glacial ice height-classified diagnostics
      real*8, dimension(:,:,:,:), allocatable :: ijhc
!@var scale_ijhc scale factor for ijhc diagnostics
      real*8, dimension(kijhc) :: scale_ijhc
!@var ia_ijhc,denom_ijhc  idacc-numbers,weights for ijhc diagnostics
      integer, dimension(kijhc) :: ia_ijhc,denom_ijhc,lgrid_ijhc
!@var sname_ijhc short names of ijhc diagnostics
      character(len=sname_strlen), dimension(kijhc) :: sname_ijhc
!@var lname_ijhc,units_ijhc descriptions/units of ijhc diagnostics
      character(len=lname_strlen), dimension(kijhc) :: lname_ijhc
      character(len=units_strlen), dimension(kijhc) :: units_ijhc
!@var cdl_ijhc consolidated metadata for ijhc output fields in cdl notation
      type(cdl_type) :: cdl_ijhc,cdl_ijhc_latlon
!@var ijhc_xxx indices for accumulations
      integer ::
     &     ijhc_frac,ijhc_fhc,ijhc_one,
     &     IJHC_SRFP,
     &     IJHC_PRECLI,IJHC_RUNLI,IJHC_EVAPLI,IJHC_F0LI,IJHC_TSLI,
     &     IJHC_SHDTLI,IJHC_EVHDT,IJHC_TRHDT,IJHC_IMPMLI,IJHC_IMPHLI

      END MODULE LANDICE_COM
! -----------------------------------------------------------
      SUBROUTINE ALLOC_LANDICE_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
#ifdef GLINT2
      use  domain_decomp_atm, only : glint2
      use glint2_modele
      use hp2hc
      use landice_com, only : fhc_approx, usedhp
#endif
      USE RESOLUTION, ONLY : IM,JM,LM
      Use LANDICE_COM, Only: NHC,FHC,
     *   SNOWLI,TLANDI, MDWNIMP,EDWNIMP,
     *   FSHGLM,FNHGLM, ELEVHP, HC_T_LAPSE_RATE
#ifdef TRACERS_WATER
      USE LANDICE_COM, ONLY : TRSNOWLI, TRLNDI, TRDWNIMP
      USE TRACER_COM, only : NTM
#ifdef TRACERS_OCEAN
      ! landice_com should really be the owner of *ACC[PDA,PDG]
      USE LANDICE, only : TRACCPDA, TRACCPDG
#endif
#endif
      use Dictionary_mod, only : sync_param, get_param
      USE LANDICE_COM, only : KIJHC,IJHC

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER

      call sync_param("HC_T_LAPSE_RATE", HC_T_LAPSE_RATE)

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

#ifdef GLINT2
      NHC = glint2_modele_nhc(glint2)
#else
      call sync_param("NHC",NHC)
#endif

      ALLOCATE( FHC(I_0H:I_1H,J_0H:J_1H,NHC),
     *          ELEVHP(I_0H:I_1H,J_0H:J_1H,NHC),
     *          SNOWLI(I_0H:I_1H,J_0H:J_1H,NHC),
     *          TLANDI(2,I_0H:I_1H,J_0H:J_1H,NHC),
     *          MDWNIMP(I_0H:I_1H,J_0H:J_1H),
     *          EDWNIMP(I_0H:I_1H,J_0H:J_1H),
     *          FSHGLM(I_0H:I_1H,J_0H:J_1H),
     *          FNHGLM(I_0H:I_1H,J_0H:J_1H),
     *          STAT=IER)
      fhc(:,:,:) = 1d0/nhc
#ifdef GLINT2
      ALLOCATE(
     *          USEDHP(I_0H:I_1H,J_0H:J_1H,NHC),
     *          FHC_APPROX(I_0H:I_1H,J_0H:J_1H,NHC))
      usedhp(:,:,:) = 0
      fhc_approx(:,:,:) = 1d0/nhc
#endif
      elevhp(:,:,:) = 0
#ifdef TRACERS_WATER
      ALLOCATE( TRSNOWLI(NTM,I_0H:I_1H,J_0H:J_1H,NHC),
     *          TRLNDI  (NTM,I_0H:I_1H,J_0H:J_1H,NHC),
     *          TRDWNIMP(NTM,I_0H:I_1H,J_0H:J_1H),
     *          STAT=IER)
#ifdef TRACERS_OCEAN
      ALLOCATE(TRACCPDA(NTM), TRACCPDG(NTM))
      TRACCPDA = 0.; TRACCPDG = 0.
#endif
#endif

      ALLOCATE(IJHC(I_0H:I_1H,J_0H:J_1H,NHC,KIJHC))

      RETURN
      END SUBROUTINE ALLOC_LANDICE_COM
! ----------------------------------------------------------

      subroutine read_landice_ic
!@sum   read_landice_ic read land ice initial conditions file.
      use model_com, only : ioread
      use domain_decomp_atm, only : grid
      use pario, only : par_open,par_close
#ifdef GLINT2
      use glint2_modele
      use domain_decomp_atm, only : glint2
      use hp2hc, only : hp_to_hc
      use landice_com, only : elevhp, fhc_approx, usedhp, fhc
      use landice_com, only : snowli, tlandi, nhc
      use pario, only : read_dist_data
      use fluxes, only : flice, flice_glint2
      use constant, only : BYGRAV
      use atm_com, only : zatmo
#endif
      use filemanager, only : file_exists
      implicit none
      integer :: fid,i

      if(file_exists('GIC')) then
        fid = par_open(grid,'GIC','read')
      elseif(file_exists('LICEIC')) then
        fid = par_open(grid,'LICEIC','read')
      else
        return
      endif

      call new_io_landice(fid,ioread)

#ifdef GLINT2
      ! Stuff above didn't read any of tlandi or snowli because
      ! Only the first (non-model) level can be expected in the GIC
      ! file.  Read the first level now.
      ! (And we'll write all HP's into the restart files).
      call read_dist_data(grid,fid,'snowli',snowli(:,:,1))
      call read_dist_data(grid,fid,'tlandi',tlandi(:,:,:,1),jdim=3)

      ! TLANDI and SNOWLI for the rest of the HP's should
      ! come out of the GLINT2 (or related) config file at
      ! some point.  For now, just copy it over and
      ! be done with it.
      do i=1,nhc
        snowli(:,:,i) = snowli(:,:,1)
        tlandi(:,:,:,i) = tlandi(:,:,:,1)
      end do
#endif
      call par_close(grid,fid)

#ifdef GLINT2
      ! FHC was just read in new_io_landice().
      ! Fix up fhc, based on glint2 API
      call glint2_modele_init_landice_com(glint2,
     &     zatmo, BYGRAV, flice_glint2, flice,
     &     usedhp, fhc, elevhp, hp_to_hc, fhc_approx,
     &     grid%i_strt_halo, grid%j_strt_halo)
#endif

      end subroutine read_landice_ic

      subroutine def_rsf_landice(fid)
!@sum  def_rsf_landice defines landice array structure in restart files
!@auth M. Kelley
!@ver  beta
      use landice_com
      use landice
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use conserv_diags
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,fhc,'fhc(dist_im,dist_jm,nhc)')
#ifdef GLINT2
      call defvar(grid,fid,usedhp,'usedhp(dist_im,dist_jm,nhc)')
      call defvar(grid,fid,fhc_approx,'fhc_approx(dist_im,dist_jm,nhc)')
#endif
      call defvar(grid,fid,elevhp,'elevhp(dist_im,dist_jm,nhc)')
      call defvar(grid,fid,snowli,'snowli(dist_im,dist_jm,nhc)')
      call defvar(grid,fid,tlandi,'tlandi(d2,dist_im,dist_jm,nhc)')
      call defvar(grid,fid,mdwnimp,'mdwnimp(dist_im,dist_jm)')
      call defvar(grid,fid,edwnimp,'edwnimp(dist_im,dist_jm)')
      call defvar(grid,fid,accpda,'accpda')
      call defvar(grid,fid,eaccpda,'eaccpda')
      call defvar(grid,fid,accpdg,'accpdg')
      call defvar(grid,fid,eaccpdg,'eaccpdg')
      call defvar(grid,fid,micbimp,'micbimp(two)')
      call defvar(grid,fid,eicbimp,'eicbimp(two)')
#ifdef TRACERS_WATER
      call defvar(grid,fid,trsnowli,
     &     'trsnowli(ntm,dist_im,dist_jm,nhc)')
      call defvar(grid,fid,trlndi,'trlndi(ntm,dist_im,dist_jm,nhc)')
      call defvar(grid,fid,trdwnimp,
     &     'trdwnimp(ntm,dist_im,dist_jm)')
#ifdef TRACERS_OCEAN
      call defvar(grid,fid,traccpda,'traccpda(ntm)')
      call defvar(grid,fid,traccpdg,'traccpdg(ntm)')
c      call defvar(grid,fid,tricbimp,'tricbimp(ntm,two)')
#endif
#endif
      call declare_conserv_diags( grid, fid, 'wlani(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'elani(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'wiceb(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'eiceb(dist_im,dist_jm)' )
      return
      end subroutine def_rsf_landice

      subroutine new_io_landice(fid,iaction)
!@sum  new_io_landice read/write landice arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      use landice_com
      use landice
      use conserv_diags
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      external conserv_MLI, conserv_MICB, conserv_HLI, conserv_HICB
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'fhc',fhc)
#ifdef GLINT2
        call write_dist_data(grid,fid,'usedhp',usedhp)
        call write_dist_data(grid,fid,'fhc_approx',fhc_approx)
#endif
        call write_dist_data(grid,fid,'elevhp',elevhp)
        call write_dist_data(grid,fid,'snowli',snowli)
        call write_dist_data(grid,fid,'tlandi',tlandi,jdim=3)
        call write_dist_data(grid,fid,'mdwnimp',mdwnimp)
        call write_dist_data(grid,fid,'edwnimp',edwnimp)
        call write_data(grid,fid,'accpda',accpda)
        call write_data(grid,fid,'eaccpda',eaccpda)
        call write_data(grid,fid,'accpdg',accpdg)
        call write_data(grid,fid,'eaccpdg',eaccpdg)
        call write_data(grid,fid,'micbimp',micbimp)
        call write_data(grid,fid,'eicbimp',eicbimp)
#ifdef TRACERS_WATER
        call write_dist_data(grid,fid,'trsnowli',trsnowli,jdim=3)
        call write_dist_data(grid,fid,'trlndi',trlndi,jdim=3)
        call write_dist_data(grid,fid,'trdwnimp',trdwnimp,jdim=3)
#ifdef TRACERS_OCEAN
        call write_data(grid,fid,'traccpda',traccpda)
        call write_data(grid,fid,'traccpdg',traccpdg)
c        call write_data(grid,fid,'tricbimp',tricbimp)
#endif
#endif
        call dump_conserv_diags( grid, fid, 'wlani', conserv_MLI )
        call dump_conserv_diags( grid, fid, 'elani', conserv_HLI )
        call dump_conserv_diags( grid, fid, 'wiceb', conserv_MICB )
        call dump_conserv_diags( grid, fid, 'eiceb', conserv_HICB )
      case (ioread)            ! input from restart file

        call read_dist_data(grid,fid,'fhc',fhc)
#ifdef GLINT2
        call read_dist_data(grid,fid,'usedhp',usedhp)
        call read_dist_data(grid,fid,'fhc_approx',fhc_approx)
#endif
        call read_dist_data(grid,fid,'elevhp',elevhp)
        call read_dist_data(grid,fid,'snowli',snowli)
        call read_dist_data(grid,fid,'tlandi',tlandi,jdim=3)
c set some defaults for quantities which may not be in the
c restart file
        mdwnimp(:,:) = 0.; edwnimp(:,:) = 0.
        accpda = 0.; eaccpda = 0.
        accpdg = 0.; eaccpdg = 0.
        MICBIMP(:) = 0  ;  EICBIMP(:) = 0
        call read_dist_data(grid,fid,'mdwnimp',mdwnimp)
        call read_dist_data(grid,fid,'edwnimp',edwnimp)
        call read_data(grid,fid,'accpda',accpda,bcast_all=.true.)
        call read_data(grid,fid,'eaccpda',eaccpda,bcast_all=.true.)
        call read_data(grid,fid,'accpdg',accpdg,bcast_all=.true.)
        call read_data(grid,fid,'eaccpdg',eaccpdg,bcast_all=.true.)
        call read_data(grid,fid,'micbimp',micbimp,bcast_all=.true.)
        call read_data(grid,fid,'eicbimp',eicbimp,bcast_all=.true.)
#ifdef TRACERS_WATER
        call read_dist_data(grid,fid,'trsnowli',trsnowli,jdim=3)
        call read_dist_data(grid,fid,'trlndi',trlndi,jdim=3)
        call read_dist_data(grid,fid,'trdwnimp',trdwnimp,jdim=3)
#ifdef TRACERS_OCEAN
        call read_data(grid,fid,'traccpda',traccpda,
     &       bcast_all=.true.)
        call read_data(grid,fid,'traccpdg',traccpdg,
     &       bcast_all=.true.)
c        call read_data(grid,fid,'tricbimp',tricbimp,bcast_all=.true.)
#endif
#endif
      end select
      return
      end subroutine new_io_landice

      subroutine def_rsf_glaacc(fid,r4_on_disk)
!@sum  def_rsf_glaacc defines accumulation array structure in restart/acc files
!@auth M. Kelley
!@ver  beta
      use landice_com, only : ijhc
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      logical :: r4_on_disk  !@var r4_on_disk if true, real*8 stored as real*4

      call defvar(grid,fid,ijhc,'ijhc(dist_im,dist_jm,nhc,kijhc)',
     &     r4_on_disk=r4_on_disk)

      return
      end subroutine def_rsf_glaacc

      subroutine new_io_glaacc(fid,iaction)
!@sum  new_io_glaacc read/write accumulation arrays from/to restart/acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use resolution, only : im,jm
      use model_com, only : ioread,iowrite,iowrite_single,idacc
      use landice_com, only : ijhc,ijhc_frac,ia_ijhc
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : hasNorthPole, hasSouthPole
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: k,l
      select case (iaction)
      case (iowrite,iowrite_single) ! output to restart or acc file
        if(iaction.eq.iowrite_single) then ! pole fills needed for acc-files
          do l=1,size(ijhc,4)
            do k=1,size(ijhc,3)
              if(hasSouthPole(grid)) then
                ijhc(2:im, 1,k,l) = ijhc(1, 1,k,l)
              endif
              if(hasNorthPole(grid)) then
                ijhc(2:im,jm,k,l) = ijhc(1,jm,k,l)
              endif
            enddo
          enddo
          do l=1,size(ijhc,4) ! mult by area fraction for masking purposes
            if(l == ijhc_frac) cycle
            ijhc(:,:,:,l) = ijhc(:,:,:,l)*
     &           (ijhc(:,:,:,ijhc_frac)/idacc(ia_ijhc(ijhc_frac)))
          enddo
        endif
        call write_dist_data(grid,fid,'ijhc',ijhc)
      case (ioread)            ! input from restart or acc file
        call read_dist_data(grid,fid,'ijhc',ijhc)
      end select
      return
      end subroutine new_io_glaacc

      subroutine def_meta_glaacc(fid)
!@sum  def_meta_glaacc defines metadata in acc files
!@auth M. Kelley
!@ver  beta
      use landice_com, only :
     &     ia_ijhc,scale_ijhc,denom_ijhc,sname_ijhc
     &     ,cdl_ijhc,cdl_ijhc_latlon
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_attr
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_attr(grid,fid,'ijhc','reduction','sum')
      call write_attr(grid,fid,'ijhc','split_dim',4)
      call defvar(grid,fid,ia_ijhc,'ia_ijhc(kijhc)')
      call defvar(grid,fid,scale_ijhc,'scale_ijhc(kijhc)')
      call defvar(grid,fid,denom_ijhc,'denom_ijhc(kijhc)')
      call defvar(grid,fid,sname_ijhc,'sname_ijhc(sname_strlen,kijhc)')
      call defvar_cdl(grid,fid,cdl_ijhc,
     &     'cdl_ijhc(cdl_strlen,kcdl_ijhc)')
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_ijhc_latlon,
     &     'cdl_ijhc_latlon(cdl_strlen,kcdl_ijhc_latlon)')
#endif

      return
      end subroutine def_meta_glaacc

      subroutine write_meta_glaacc(fid)
!@sum  write_meta_glaacc writes metadata to acc files
!@auth M. Kelley
!@ver  beta
      use landice_com, only :
     &     ia_ijhc,scale_ijhc,denom_ijhc,sname_ijhc
     &     ,cdl_ijhc,cdl_ijhc_latlon
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_data
      use cdl_mod, only : write_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_data(grid,fid,'ia_ijhc',ia_ijhc)
      call write_data(grid,fid,'scale_ijhc',scale_ijhc)
      call write_data(grid,fid,'denom_ijhc',denom_ijhc)
      call write_data(grid,fid,'sname_ijhc',sname_ijhc)
      call write_cdl(grid,fid,'cdl_ijhc',cdl_ijhc)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_ijhc_latlon',cdl_ijhc_latlon)
#endif

      return
      end subroutine write_meta_glaacc
