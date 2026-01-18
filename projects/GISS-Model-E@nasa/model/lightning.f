#include "rundeck_opts.h"
  
      module lightning
!@sum  lightning variables for lightning parameterization
!@auth Colin Price (modification by Greg Faluvegi and Lee Murray)
      USE RESOLUTION, ONLY : LM

      implicit none
      save

      ! Set default parameterization for lightning flash rate
      INTEGER,             PUBLIC :: lightning_param = 1
      ! 1 = Cloud Top Height Scheme (default) [Price and Rind, 1992, 1993, 1994]
      ! 2 = Convective Mass Flux Scheme [Allen et al., 2001]
      ! 3 = Convectiev Precipitation Scheme [Allen and Pickering, 2002]

#ifdef TRACERS_SPECIAL_Shindell
      REAL*8, ALLOCATABLE, PUBLIC :: CLDTOPL(:,:)   ! Cloud-top level
      REAL*8, ALLOCATABLE, PUBLIC :: ENOx_lgt(:,:)  ! Col total src [kg(N)/s]
      INTEGER, PARAMETER,  PUBLIC :: NLTYPE  = 4    ! Lightning profile types
      INTEGER, PARAMETER,  PUBLIC :: NNLIGHT = 3200 ! Number of pts in profile
      REAL*8, ALLOCATABLE, PUBLIC :: LNOX_CDF(:,:)  ! Ott et al. [2010] CDFs

      ! NOx yields per flash type (moles N/flash)
      ! Initialize with value that reproduces methane lifetime in present-day
      ! and minimizes tropospheric ozone bias (160 molN/flash)
      REAL*8,              PUBLIC  :: FLASH_YIELD_MIDLAT = 160d0 ! moles N/fl
      REAL*8,              PUBLIC  :: FLASH_YIELD_TROPIC = 160d0 ! moles N/fl

#endif

      ! SUBDD Instantaneous Diagnostic Arrays
      REAL*8, ALLOCATABLE, PUBLIC :: FLASH_DENS(:,:)! Flash density [fl/m2/s]
      REAL*8, ALLOCATABLE, PUBLIC :: CG_DENS(:,:)   ! Cloud-to-ground density

#ifdef AUTOTUNE_LIGHTNING
      ! If AUTOTUNE_LIGHTNING is defined in the rundeck, then the tuning
      ! parameters (which are dependent on horizontal resolution, and very
      ! sensitive to any changes in the convective cloud code) are determined
      ! online by comparison with the long-term climatology from the LIS/OTD
      ! satellites. Run the model at least a couple years, and discard from
      ! analysis any years before the tuning parameters printed in the PRT
      ! file stabalize. This saves the hassle from having to run the model to
      ! generate a climatology to compare to the satellites, and then re-run
      ! with tuning parameters prescribed. This option is only appropriate
      ! for present-day simulations, but values calculated from present-day 
      ! runs may then be used in future and paleo atmospheres.
      ! Contact Lee T. Murray with any questions.

      ! Historical Flash Rates for Determining Scaling Factor
      INTEGER, PUBLIC, PARAMETER  :: NHISTLI = 87600
      REAL*8, ALLOCATABLE, PUBLIC :: FLASH_UNC(:,:) ! Unconstrained SUBDD dens
      REAL*8,  PUBLIC             :: LAND_FR_UNC(NHISTLI) = 0d0 ! [flash/s]
      REAL*8,  PUBLIC             ::  SEA_FR_UNC(NHISTLI) = 0d0 ! [flash/s]
      REAL*8,  PUBLIC             ::      CNT_FR(NHISTLI) = 0d0 ! [count] 

      ! OTD/LIS Climatological Values [flash/s]
      REAL*8,  PUBLIC, PARAMETER  :: LAND_FR_LIS = 37.81122d0
      REAL*8,  PUBLIC, PARAMETER  ::  SEA_FR_LIS = 12.18967d0
#endif

      REAL*8,  PARAMETER, PUBLIC  :: T_NEG_BOT     = 273.0d0  !   0 C 
      REAL*8,  PARAMETER, PUBLIC  :: T_NEG_CTR     = 258.0d0  ! -15 C
      REAL*8,  PARAMETER, PUBLIC  :: T_NEG_TOP     = 233.0d0  ! -40 C

      ! Scaling parameters
      ! These make the present-day climatological mean lightning flash rate 
      ! match the LIS/OTD satellite observations (~50 flashes/s).
      ! These are determiend for F40 TCadi, reset in rundeck for others
#if defined( LIGHTNING_MFLUX )
      REAL*8,             PUBLIC :: tune_lt_land = 1.65d0
      REAL*8,             PUBLIC :: tune_lt_sea  = 0.29d0
#elif defined( LIGHTNING_PRECON )
      REAL*8,             PUBLIC :: tune_lt_land = 9.61d0
      REAL*8,             PUBLIC :: tune_lt_sea  = 0.60d0
#else
      REAL*8,             PUBLIC :: tune_lt_land = 2.31d0
      REAL*8,             PUBLIC :: tune_lt_sea  = 5.52d0
#endif

      ! Set pertubations to flash rate in rundeck, ratio
      REAL*8,             PUBLIC  :: FLASH_PERTURB = 1d0

      ! model level one below the one closest to nominal mid-level
      ! pressure 440 hPa:
      INTEGER :: L440mbM1=0

      end module lightning

      subroutine alloc_lightning(grid)

!@SUM  alllocate lightning arrays for current grid
!@auth L. Murray
!@ver  1.0

      USE DOMAIN_DECOMP_ATM, ONLY : dist_grid, getDomainBounds
      USE RESOLUTION,        ONLY : LM
      USE ATM_COM,           ONLY : PMIDL00
      USE LIGHTNING,         ONLY : CG_DENS
      USE LIGHTNING,         ONLY : FLASH_DENS
      USE LIGHTNING,         ONLY : L440mbM1
#ifdef TRACERS_SPECIAL_Shindell
      USE FILEMANAGER,       ONLY : OPENUNIT, CLOSEUNIT
      USE LIGHTNING,         ONLY : ENOx_lgt, CLDTOPL
      USE LIGHTNING,         ONLY : LNOx_CDF
      USE LIGHTNING,         ONLY : NNLIGHT
      USE LIGHTNING,         ONLY : NLTYPE
#endif
#ifdef AUTOTUNE_LIGHTNING
      USE LIGHTNING,         ONLY : FLASH_UNC
#endif
      implicit none

      type (dist_grid), intent(in) :: grid

      INTEGER             :: J_1H, J_0H, I_0H, I_1H
      INTEGER             :: AS, III, IOS, JJJ, IU_FILE, L
      REAL*8              :: Y0, Y1
      CHARACTER(LEN=255)  :: FILENAME

      !======================================
      ! ALLOC_LIGHTNING begins here
      !======================================

      ! Get local grid info
      call getDomainBounds(grid, 
     &               I_STRT_HALO=I_0H, I_STOP_HALO=I_1H,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      !------------------------------------
      ! Allocate arrays for lightning
      !------------------------------------

      allocate( FLASH_DENS(I_0H:I_1H,J_0H:J_1H) )
      FLASH_DENS = 0d0

      allocate(    CG_DENS(I_0H:I_1H,J_0H:J_1H) )
      CG_DENS = 0d0

#ifdef AUTOTUNE_LIGHTNING
      allocate( FLASH_UNC(I_0H:I_1H,J_0H:J_1H) )
      FLASH_UNC = 0d0
#endif
      !---------------------------------------------
      ! Define nominal model level one below the one
      ! closest to 440 hPa
      !---------------------------------------------
      L440mbM1=minloc(abs(PMIDL00-440.d0),1)-1
      if(L440mbM1.lt.1.or.L440mbM1.gt.LM)call stop_model('L440mbM1',255)

#ifdef TRACERS_SPECIAL_Shindell
 
      !------------------------------------
      ! Allocate arrays for NOx emissions
      !------------------------------------

      allocate(      ENOx_lgt(I_0H:I_1H,J_0H:J_1H) )
      allocate(      CLDTOPL(I_0H:I_1H,J_0H:J_1H)  )
      allocate(      LNOx_CDF( NNLIGHT, NLTYPE )   )

      ENOx_lgt   = 0d0
      CLDTOPL    = 0d0
      LNOx_CDF   = 0d0

      !-----------------------------------------------------------------
      ! Read lightning CDF from Ott et al [JGR, 2010]. (ltm, 1/25/11)
      !-----------------------------------------------------------------

      ! Open file containing lightning PDF data
      CALL OPENUNIT( 'LNOxCDF', IU_FILE, .false., .true. )

      ! Read 12 header lines
      DO III = 1, 12
         READ( IU_FILE, '(a)', IOSTAT=IOS ) 
         !IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:2' )
      ENDDO

      ! Read NNLIGHT types of lightning profiles
      DO III = 1, NNLIGHT
         READ( IU_FILE,*,IOSTAT=IOS) (LNOx_CDF(III,JJJ),JJJ=1,NLTYPE)
      ENDDO

      ! Close file
      CALL CLOSEUNIT( IU_FILE )

#endif

      return
      end subroutine alloc_lightning
     
      subroutine calc_lightning(i,j,lmax,lfrz,mflux,precon)
!@sum calc_lightning calculates lightning flash amount and cloud-
!@+   to-ground amount, based on cloud top height.

      USE LIGHTNING, ONLY : lightning_param
      USE LIGHTNING, ONLY : TUNE_LT_LAND, TUNE_LT_SEA, FLASH_PERTURB
      USE LIGHTNING, ONLY : CG_DENS, FLASH_DENS, T_NEG_CTR
#ifdef TRACERS_SPECIAL_Shindell
      USE LIGHTNING, ONLY : FLASH_YIELD_MIDLAT, FLASH_YIELD_TROPIC
      USE LIGHTNING, ONLY : CLDTOPL, ENOx_lgt
#endif
#ifdef AUTOTUNE_LIGHTNING
      use DOMAIN_DECOMP_1D, only: am_i_root
      USE LIGHTNING, ONLY : FLASH_UNC
#endif
      use atm_com,           only : zatmo, gz
      use model_com,         only : DTsrc
      use geom,              only : lat2d_dg,axyp,byaxyp
      use constant,          only : bygrav
      use Dictionary_mod,    only : sync_param
      use TimeConstants_mod, only : SECONDS_PER_DAY
      use fluxes,            only : focean, atmsrf ! contains tsavg
      use diag_com,          only : ij_CtoG,ij_flash,aij=>aij_loc
 
      implicit none

!@var lmax highest layer of current convective cloud event
!@var lfrz freezing level
!@var LTOP local copy of LMAX (alterable)
!@var htcon height of convection?
!@var htfrz height of freezing?
!@var flash lightning flashes per minute
!@var th,th2,th3,th4 thickness of cold sector, squared, cubed, etc.
!@var cg fraction of lightning that is cloud-to-ground
!@var zlt is a ratio of IC/CG see Price and Rind (1993) for more
      integer, intent(in) :: lmax,lfrz,i,j
      real*8,  intent(in) :: mflux, precon
      integer:: LTOP, L
      real*8 :: htcon,htfrz,flash,th,th2,th3,th4,zlt,cg,area_ref
      real*8 :: htcon2,flashun,pr,mf,cgun


      ! Initialize
      flash   = 0d0
      cg      = 0d0
      flashun = 0d0
      cgun    = 0d0

      call sync_param("lightning_param",lightning_param)
#ifdef TRACERS_SPECIAL_Shindell
      call sync_param("FLASH_YIELD_MIDLAT",FLASH_YIELD_MIDLAT)
      call sync_param("FLASH_YIELD_TROPIC",FLASH_YIELD_TROPIC)
#ifndef AUTOTUNE_LIGHTNING
      call sync_param("tune_lt_land",tune_lt_land)
      call sync_param("tune_lt_sea",tune_lt_sea)
#else 
      if ( am_I_root() )  THEN
         WRITE(6,*) FLASH_YIELD_MIDLAT, FLASH_YIELD_TROPIC
         WRITE(6,*) tune_lt_land, tune_lt_sea
      end if
#endif
#endif

      ! Level of highest convective activity
      ltop=lmax
      
      ! If surface air temperature is already colder than the negative
      ! charge layer of the cloud (-15C [Houze,1993]), then assume no lightning.

      if ( atmsrf%tsavg(i,j) .gt. T_NEG_CTR ) then

         !===============================================
         ! (1) Calculate the CG/IC ratio
         !===============================================

         ! Geopotential (gz) is relative to mean sea level,
         ! so subtract surface elevation (zatmo/g).
         htcon=max((gz(i,j,ltop)-zatmo(i,j))*bygrav,0d0)*1d-3
         htfrz=max((gz(i,j,lfrz)-zatmo(i,j))*bygrav,0d0)*1d-3

         ! Calculate the fraction of total lightning that is cloud-to-ground
         ! using Price and Rind [GRL, 1993]. This uses the thickness of the
         ! cold sector of the cloud (Hmax-Hzero) as the determining
         ! parameter.  The algorithm is only valid for thicknesses greater
         ! than 5.5km and less than 14km.
         
         th=(htcon-htfrz)
         th=min(max(th,5.5d0),14.d0)
         th2=th*th
         th3=th2*th
         th4=th3*th
         zlt = 0.021d0*th4 - 0.648d0*th3 + 7.493d0*th2 
     &        - 36.544d0*th + 63.088d0
         

         !===============================================
         ! (2) Calculate the flash rate (flashes/s)
         !===============================================

         ! Which parameterization are we using?
         select case ( lightning_param ) 

         case ( 2 ) ! Convective mass flux

         !===============================================
         ! (2a) Convective Mass Flux Scheme
         !===============================================
         ! Allen and Pickering (2002) give the following 
         ! parameterizations for lightning flash rates as a function 
         ! of upward cloud mass flux [kg m^-2 min^-1] at 0.44 sigma (440mb)
         !
         !    LF_CG = [delta x][delta y] *
         !            ( a + b*M + c*M^2 + d*M^3 + e*M^4 ) / A
         !
         !    For: 0 < M < 9.6 [kg/m2/min]
         !
         ! Where:
         !    (1) LF_CG is the CG flash rate [flashes/min)] within the
         !         2.0 x 2.5 grid box 
         !    (2) a, b, c, d, e are coefficients, listed below
         !    (3) [delta x][delta y] is the area of the grid box
         !    (4) A is the area of a 2.0 x 2.5 box centered at 30N
         !    (5) M is the upward cloud mass flux at 0.44 sigma
         !        
         ! Since the polynomial experiences an inflection at 
         ! M ~= 9.6 [kg/m2/min], points greater than this are 
         ! set to 9.6 [kg/m2/min].
         ! 
         ! The polynomial coefficients are:
         !    a=-2.34e-2, b=3.08e-1, c=-7.19e-1, d=5.23e-1, e=-3.71e-2
         !==============================================================

            ! Convert from hPa/m2 to kg/m2/min
            mf = mflux * ( 100 * bygrav ) * ( 60d0 / DTsrc )
         
            ! The first equation is a unified equation for all surface types
            ! The second two equations follow the general land/sea dichotomy

            ! Limit mass flux range from 0.9-9.6 [kg/m2/min]
            mf=min(max(mf,0.096d0),9.6d0)

            ! Calculate cloud-to-ground flash rate
            cgun = -2.34d-02 + mf * (  3.08d-01 +
     &                         mf * ( -7.19d-01 +
     &                         mf * (  5.23d-01 +
     &                         mf * ( -3.71d-02 ) ) ) )

            ! Disallow negative values that the polynomial can produce 
            cgun = max( cgun, 0d0 )

            ! Normalize by the area of a grid box at 30N
            area_ref = 5.3528d10 ! avg 2x2.5 box at 30N
            cgun = cgun * ( axyp(i,j) / area_ref )
         
            ! Convert from flashes/min to flashes/s
            cgun = cgun / 60d0

            ! Calculate total flash rate from cloud-to-ground rate
            flashun = cgun * (1.+zlt)

            ! Apply tuning parameter
            flash = (      focean(i,j)  * tune_lt_sea  + 
     &                (1d0-focean(i,j)) * tune_lt_land   ) * flashun
            cg    = flash / (1.+zlt)

         case ( 3 ) ! Convective Precipitation Scheme

         !===============================================
         ! (2b) Convective Precipitation
         !===============================================
         ! Allen and Pickering (2002) give the following parameterizations 
         ! for CG lightning flash rates as a function of surface 
         ! convective precipitation [mm d^-1]
         !
         !    LF_CG = [delta x][delta y] *
         !            ( a + b*P + c*P^2 + d*P^3 + e*P^4 ) /A,
         !                                           
         !    For: 7 < P < 90 [mm/day]
         !       
         ! Where:
         !    (1) LF_CG = CG flash rate (flashes/min) w/in the grid box
         !    (2) a, b, c, d, e are coefficients, listed below
         !    (3) [delta x][delta y] is the area of the grid box
         !    (4) A is the area of a grid box centered at 30N
         !    (5) P is the PRECON rate [mm/day] during the 6-hour period.
         !
         ! The polynomial coefficients for land boxes are:
         !   a=3.75e-02, b=-4.76e-02, c=5.41e-03, d=3.21e-04, e=-2.93e-06
         !
         ! and the polynomial coefficients for water boxes are:
         !   a=5.23e-02, b=-4.80e-02, c=5.45e-03, d=3.68e-05, e=-2.42e-07
         !
         ! Both polynomials give slightly negative values for small precip
         ! rates. In addition, the land polynomial has an inflection point 
         ! at ~90 [mm/day]. Therefore flash rates are set to 0 for precip
         ! rates of less than 7 [mm/day] and to the value at 90 [mm/day]
         ! precipitation rates exceeding 90 [mm/day].
         !=================================================================

            ! Convert from hPa/m2 to mm/day, Note: 1kg(H2O)/m2 = 1mm H2O
            pr = precon * ( 100 * bygrav ) * ( SECONDS_PER_DAY / DTsrc )

            ! Limit pr range to 5.74-90 mm/day
            pr=min(max(pr,5.74d0),90d0) 
         
            cgun = focean(i,j) *
     &              5.23d-02 + pr * ( -4.80d-02 +
     &                         pr * (  5.45d-03 +
     &                         pr * (  3.68d-05 +
     &                         pr * ( -2.42d-07 ) ) ) ) + 
     &      (1d0 - focean(i,j)) *
     &              3.75d-02 + pr * ( -4.76d-02 +
     &                         pr * (  5.41d-03 +
     &                         pr * (  3.21d-04 +
     &                         pr * ( -2.93d-06 ) ) ) )

            ! Disallow negative values
            cgun = max( cgun, 0d0 )

            ! Normalize by the area of a grid box at 30N
            area_ref = 5.3528d10 ! avg 2x2.5 box at 30N
            cgun = cgun * ( axyp(i,j) / area_ref )

            ! Convert from flashes/min to flashes/s
            cgun = cgun / 60d0

            ! Calculate total flash rate from cloud-to-ground rate
            flashun = cgun * (1.+zlt)
         
            ! Apply tuning parameter
            flash = (      focean(i,j)  * tune_lt_sea  + 
     &                (1d0-focean(i,j)) * tune_lt_land   ) * flashun
            cg    = flash / (1.+zlt)

         case default ! Cloud-top height scheme (default)

         !===============================================
         ! (2c) Cloud-Top Height scheme
         !===============================================
         ! Refs: Price and Rind [1992, 1993, 1994]

            flashun = (focean(i,j)) * 6.40d-4*(htcon**1.73d0) + ! ocean
     &            (1d0-focean(i,j)) * 3.44d-5*(htcon**4.92d0)   ! land

#ifdef CUBED_SPHERE
            ! rescale flash rate by gridbox area
            area_ref = 6d10     ! avg 2x2.5 tropical area for reference
            flashun  = flashun * axyp(i,j)/area_ref
#endif

            ! Apply tuning parameter
            flash = (      focean(i,j)  * tune_lt_sea  + 
     &                (1d0-focean(i,j)) * tune_lt_land   ) * flashun

            ! Convert from flashes/min to flashes/s
            flashun = flashun / 60d0
            flash   = flash   / 60d0

            ! Calculate cloud-to-ground flash rate from total flash rate
            cg    = flash / (1.+zlt)
            
         end select

! End of different flash rate parameterization options

         ! Apply any perturbations to flash rate
         flash = flash * FLASH_PERTURB
         cg    = cg    * FLASH_PERTURB

      end if

      !===============================================
      ! (3) Archive Diagnostics
      !===============================================

      ! Accumulate lightning flashes as flashes/m2 (convert from flashes/s)
      aij(i,j,ij_flash) = aij(i,j,ij_flash) + flash*DTsrc*byaxyp(i,j)
      aij(i,j,ij_CtoG)  = aij(i,j,ij_CtoG)  +    cg*DTsrc*byaxyp(i,j)

      ! Also save for SUBDD instantaneous output
      FLASH_DENS(i,j)   =   flash*byaxyp(i,j)  ! flashes/s -> flashes/m2/s
      CG_DENS(i,j)      =      cg*byaxyp(i,j)  ! flashes/s -> flashes/m2/s
#ifdef AUTOTUNE_LIGHTNING
      FLASH_UNC(i,j)    =             flashun  ! flashes/s 
#endif

#ifdef TRACERS_SPECIAL_Shindell

      !===============================================
      ! (4) Generate NOx yield from lightning_paraming
      !===============================================
      
      ! Assume higher yields in the extratropics than the tropics following
      ! in situ and remote observational evidence [Schumann & Huntrieser, 2007]

      ! Note: GISS model used to assume IC flashes generated 10% the NOx of a
      ! CG flash. Remove that assumption following Ott et al. [2010]

      if ( abs(lat2d_dg(i,j)) > 23 ) then
         ! fl/s*mol/fl*kg(N)/mol -> kg(N)/s
         ENOx_lgt(i,j) = flash * FLASH_YIELD_MIDLAT * 14d-3 
      else
         ENOx_lgt(i,j) = flash * FLASH_YIELD_TROPIC * 14d-3
      endif

      ! Save cloud-top level for distributing NOx vertically
      CLDTOPL(i,j) = LTOP

#endif

      end subroutine calc_lightning

#ifdef TRACERS_SPECIAL_Shindell
! !DESCRIPTION: Subroutine LIGHTDIST uses a CDF to distribution the total 
!  column lightning NOx into the model vertical layers.
!
! !INTERFACE:
!
      SUBROUTINE LIGHTDIST( I, J, LTOP, dZ, TOTAL, VERTPROF )
!
! !USES:
!
      use model_com,           only : modelEclock
      use resolution,          only : LM
      use fluxes,              only : focean
      use lightning,           only : lnox_cdf, nnlight
      use geom,                only : lat2d_dg 
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I            ! Longitude index
      INTEGER,        INTENT(IN)  :: J            ! Latitude index 
      INTEGER,        INTENT(IN)  :: LTOP         ! Level of conv cloud top
      REAL*8,         INTENT(IN)  :: dZ(LM)       ! Thickness [m] < LTOP; 
                                                  !           0d0 above LTOP 
      REAL*8,         INTENT(IN)  :: TOTAL        ! Column Total LNOx molec 
!
! !OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(OUT) :: VERTPROF(LM) ! Vertical profile of LNOx
!
! !LOCAL VARIABLES:
!
      INTEGER            :: M, MTYPE, L, III, IOS, IUNIT, JJJ, JMON
      REAL*8             :: ZH, CTH
      REAL*8             :: FRAC(LM)
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! LIGHTDIST begins here!
      !=================================================================

      ! Initialize 
      FRAC     = 0d0
      MTYPE    = 0
      VERTPROF = 0d0

      ! Convective cloud top height
      CTH      = SUM(dZ)

      ! Current month
      JMON     = modelEclock%getMonth()

      !=================================================================
      ! Test whether location (I,J) is continental or marine
      !
      ! Depending on the combination of land/water and latitude, 
      ! assign a flag describing the lightning profile to apply:
      !
      !   MTYPE = 1: ocean lightning
      !   MTYPE = 2: tropical continental lightning
      !   MTYPE = 3: midlatitude continental lightning 
      !   MTYPE = 4: subtropical lightning
      !             
      ! (ltm, bmy, 1/25/11)
      !=================================================================

      ! Assign profile kind to grid box, following Allen et al. [JGR, 2010] 
      ! (ltm, 1/25/11)
      SELECT CASE ( JMON )
         
      ! Southern Hemisphere Summer
      CASE ( 1,2,3,12 )
         
         IF ( ABS( lat2d_dg(i,j) ) .le. 15 ) THEN
            IF ( focean(i,j) < 0.5 ) THEN
               MTYPE = 2        ! Tropical continental
            ELSE
               MTYPE = 1        ! Tropical marine
            ENDIF
         ELSE IF ( ( lat2d_dg(i,j) .gt. 15. ) .and. 
     &             ( lat2d_dg(i,j) .le. 30. )        ) THEN
            MTYPE = 4           ! N. Subtropics
         ELSE IF ( ( lat2d_dg(i,j) .ge. -40. ) .and. 
     &             ( lat2d_dg(i,j) .lt. -15. )       ) THEN
            MTYPE = 4           ! S. Subtropics
         ELSE
            MTYPE = 3           ! Midlatitude
         ENDIF
         
      ! Equinox months
      CASE ( 4,5,10,11 )
         
         IF ( ABS( lat2d_dg(i,j) ) .le. 15 ) THEN
            IF ( focean(i,j) < 0.5 ) THEN
               MTYPE = 2        ! Tropical continental
            ELSE
               MTYPE = 1        ! Tropical marine
            ENDIF
         ELSE IF ( ABS(lat2d_dg(i,j) ) .le. 30 ) THEN
            MTYPE = 4           ! Subtropics
         ELSE
            MTYPE = 3           ! Midlatitude
         ENDIF
         
      ! Northern Hemisphere Summer
      CASE ( 6,7,8,9 )
         
         IF ( ABS( lat2d_dg(i,j) ) .le. 15 ) THEN
            IF ( focean(i,j) < 0.5 ) THEN
               MTYPE = 2        ! Tropical continental
            ELSE
               MTYPE = 1        ! Tropical marine
            ENDIF
         ELSE IF ( ( lat2d_dg(i,j) .gt.  15. ) .and. 
     &             ( lat2d_dg(i,j) .le.  40. )       ) THEN
            MTYPE = 4           ! N. Subtropics
         ELSE IF ( ( lat2d_dg(i,j) .ge. -30. ) .and. 
     &             ( lat2d_dg(i,j) .lt. -15. )       ) THEN
            MTYPE = 4           ! S. Subtropics
         ELSE
            MTYPE = 3           ! Midlatitude
         ENDIF
         
      END SELECT

      !=================================================================
      ! Use the CDF for this type of lightning to partition the total
      ! column lightning into the model layers
      !=================================================================

      ! Compute the height [km] at the top of each vertical level.
      ! Look up the cumulative fraction of NOx for each vertical level

      ! Reinitialize height
      ZH = 0.0

      DO L = 1, LTOP

         ! Height above ground of grid box top [m]
         ZH = ZH + dZ(L)

         ! Fraction of lightning released beneath that height [unitless]
         FRAC(L) = LNOX_CDF( NINT((ZH/CTH)*REAL(NNLIGHT)), MTYPE )*0.01

      ENDDO

      ! Convert from cumulative fraction to fraction for each level
      DO L = LTOP, 2, - 1
         FRAC(L) = FRAC(L) - FRAC(L-1)
      ENDDO 
  
      ! Partition lightning NOx by layer into VERTPROF
      DO L = 1, LTOP
         VERTPROF(L) = ( FRAC(L) * TOTAL )
      ENDDO

      END SUBROUTINE LIGHTDIST

      subroutine get_lightning_NOx
!@sum  get_lightning_NOx to define the 3D source of NOx from lightning
!@auth Colin Price / Greg Faluvegi / Lee Murray
 
      use geom, only                : lat2d_dg,byaxyp
      use fluxes, only              : tr3Dsource
      use tracer_com, only          : n_NOx,nOther
      use lightning, only           : ENOx_lgt, cldtopl
      use constant, only            : bygrav
      use atm_com, only             : zatmo, gz, ltropo
      use resolution, only          : LM
      use domain_decomp_atm, only   : GRID, getDomainBounds
#ifdef ACCMIP_LIKE_DIAGS
      use trdiag_com, only          : taijls=>taijls_loc,ijlt_NOxLgt
#endif
 
      implicit none
 
      real*8, dimension(LM) :: dZ
      real*8, dimension(LM) :: SRCLIGHT 
      integer               :: i, j, l, ltop, levtrop

      integer               :: J_1, J_0, J_0H, J_1H, I_0, I_1

      I_0 = grid%I_STRT; J_0 = grid%J_STRT
      I_1 = grid%I_STOP; J_1 = grid%J_STOP
 
      do j=J_0,J_1
      do i=I_0,I_1

         ! Reset arrays
         srclight(:)                    = 0d0
         tr3Dsource(i,j,:,nOther,n_NOx) = 0d0

         ! Skip columns with no lightning emissions
         if ( ENOx_lgt(i,j) .eq. 0d0 ) cycle

         ! Set some variables
         ltop    = cldtopl(i,j) ! Local cloud top level
         levtrop = ltropo(i,j)  ! Local tropopause level 

         ! Get vertical thickness of each column [m]
         dZ(:) = 0d0
         if ( ltop < LM ) then
            do l=1,ltop
               dZ(l) = (0.5*(gz(i,j,l)+gz(i,j,l+1))-zatmo(i,j))*bygrav
               ! Until I can figure out how to get real temperature
               !dZ(l) = rgas/grav*TX(i,j,L)*log(pedn(l,i,j)/pedn(L+1,i,j))
            enddo
         end if
         !if ( sum(dZ) == 0 ) cycle ! This may happen on the first call.

         ! Distribute ENOx_lgt column emissions vertically into srclight
         !   using cumulative distribution functions from Ott et al. [2010]
         call lightdist( i, j, ltop, dZ, ENOx_lgt(i,j), srclight )

         ! Save tracer 3D source. Note: This puts any stratospheric LNOx into 
         ! the UT. Check to see how systematic this is, and whether or not we 
         ! should put it in the LS.
         do L=1,levtrop
            tr3Dsource(i,j,      L,nOther,n_NOx) = srclight(L)
         enddo 
         do L=levtrop+1,LM
            tr3Dsource(i,j,levtrop,nOther,n_NOx) = 
     &      tr3Dsource(i,j,levtrop,nOther,n_NOx) + srclight(L)
         enddo  
#ifdef ACCMIP_LIKE_DIAGS
         do L=1,LM
           taijls(i,j,L,ijlt_NOxLgt)=taijls(i,j,L,ijlt_NOxLgt) +
     &     tr3Dsource(i,j,L,nOther,n_NOx)*byaxyp(i,j)
         enddo
#endif

      enddo  ! I
      enddo  ! J
         
      end subroutine get_lightning_NOx
#endif

#ifdef AUTOTUNE_LIGHTNING

      subroutine def_rsf_lightning(fid)
!@sum  def_rsf_lightning defines lightning array structure in restart files
!@auth L. Murray
      use lightning, only : LAND_FR_UNC, SEA_FR_UNC, CNT_FR, NHISTLI
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid !@var fid file id

      call defvar(grid,fid,land_fr_unc,'land_fr_unc(nhistli)')
      call defvar(grid,fid,sea_fr_unc,'sea_fr_unc(nhistli)')
      call defvar(grid,fid,cnt_fr,'cnt_fr(nhistli)')

      return
      end subroutine def_rsf_lightning

      subroutine new_io_lightning(fid,iaction)
!@sum  new_io_clouds read/write lightning arrays from/to restart files
!@auth L. Murray
      use lightning, only :  LAND_FR_UNC, SEA_FR_UNC, CNT_FR
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_data,read_data
      implicit none
      integer fid     !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)  ! output to restart file
         call write_data(grid, fid, 'land_fr_unc', land_fr_unc )
         call write_data(grid, fid, 'sea_fr_unc',   sea_fr_unc )
         call write_data(grid, fid, 'cnt_fr',           cnt_fr )
      case (ioread)   ! input from restart file
         call read_data(grid, fid, 'land_fr_unc', land_fr_unc )
         call read_data(grid, fid, 'sea_fr_unc',   sea_fr_unc )
         call read_data(grid, fid, 'cnt_fr',           cnt_fr )
      end select
      return
      end subroutine new_io_lightning

#endif
