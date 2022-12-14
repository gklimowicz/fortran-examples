#include "rundeck_opts.h"

      module AerParam_mod
!@sum This module reads, time-interpolates, and stores fields needed
!@+   by the radiation code in the prescribed-aerosol configuration
!@+   of modelE.  Subroutine updateAerosol2 provides the fields
!@+   used to calculate the direct radiative effects of aerosols.
!@+   Subroutine dCDNC_est provides a parameterized estimate of
!@+   changes of lower-tropospheric CDNC relative to 1850, as input
!@+   to prescriptions of aerosol indirect effects on cloud cover
!@+   and optical depth.
!@+   Dust aerosols are currently ingested via a separate module.
!@auth D. Koch, R. Ruedy
!@auth M. Kelley added comments, reprogrammed for netcdf input

      use timestream_mod, only : timestream
      implicit none
      save
      private

      public :: dCDNC_est
      public :: updateAerosol
      public :: updateAerosol2
      public :: DRYM2G, aermix
      public :: lma

      real*8, allocatable :: anssdd(:,:)
      real*8, allocatable :: mdpi(:,:,:)
      real*8, allocatable :: mdcur(:,:,:)
      real*8, allocatable :: md1850(:,:,:,:)

      character(len=3), dimension(6), parameter :: aernames=(/
!**** Sulfate
     &     'SUL',
!**** Sea Salt
     &     'SSA',
!**** Nitrate
     &     'NIT',
!**** Organic Carbon
     &     'OCA',
!**** Black Carbon from Fossil and bio fuel
     &     'BCA',
!**** Black Carbon from Biomass burning
     &     'BCB'
     &     /)

      real*8 :: DRYM2G(8) =
     &     (/4.667, 0.866, 4.448, 5.017, 9.000, 9.000, 1.000,1.000/)

C     Layer  1    2    3    4    5    6    7    8    9
      INTEGER :: La720=3 ! top low cloud level (aerosol-grid).
                         ! =3 for orig 9-level model
      REAL*8 , PARAMETER ::
     &     Za720=2635.                ! depth of low cloud region (m)
     &    ,byz_cm3 = 1.d-6 / Za720    ! 1d-6/depth in m (+conversion /m3 -> /cm3)
     &    ,byz_gcm3 = 1.d-3 * byz_cm3 ! g vs kg

      REAL*8, dimension(13) :: AERMIX=(/
C      Pre-Industrial+Natural 1850 Level  Industrial Process  BioMBurn
C      ---------------------------------  ------------------  --------
C       1    2    3    4    5    6    7    8    9   10   11   12   13
C      SNP  SBP  SSP  ANP  ONP  OBP  BBP  SUI  ANI  OCI  BCI  OCB  BCB
     + 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 1.9, 1.0, 1.0, 2.5, 1.9, 2.5, 1.9/)

      integer :: ima, jma, lma

!@var A6streams interface for reading and time-interpolating AERO files
!@+   See usage notes in timestream_mod
      type(timestream), dimension(6) :: A6streams

#ifdef OLD_BCdalbsn
!@var BCdepstream interface for reading and time-interpolating BC_dep file
      type(timestream), public :: BCdepstream
!@var depoBC,depoBC_1990 prescribed black carbon deposition (curr,1990)
!@+   for parameterization of the BC effect on snow albedo
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: depoBC,depoBC_1990
#else
!@var BCdalbsnstream interface for reading and time-interpolating BC_dalbsn file
      type(timestream), public :: BCdalbsnstream
!@var BCdalbsn prescribed delta-albedo of snow (units = 1) on land/seaice from BC effects
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: BCdalbsn
#endif

      contains

      subroutine dCDNC_EST(i,j,pland, dCDNC) !, table)
!@sum  finds change in cloud droplet number concentration since 1850
!@auth R. Ruedy
!@ver  1.0
      USE CONSTANT, only : pi
      implicit none
      integer, intent(in)  :: i,j ! grid indices
      real*8 , intent(in)  :: pland ! land fraction
      real*8 , intent(out) :: dCDNC ! CDNC(cur)-CDNC(1850)

      real*8, parameter, dimension(5) ::
C                TROPOSPHERIC AEROSOL PARAMETERS
C                  SO4     NO3    OCX    BCB   BCI
     &  f_act=(/ 1.0d0,  1.0d0, 0.8d0, 0.6d0, .8d0/), ! soluble fraction
     &  dens =(/1769d0, 1700d0,  1.d3,  1.d3, 1.d3/)  ! density

      real*8, parameter, dimension(2) ::
C                    Ocean         Land      ! r**3: r=.085,.052 microns
     &  radto3 =(/ 614.125d-24, 140.608d-24/),  ! used for SO4,NO3,OC,BC
     &  scl    =(/     162d0,       298d0/),  ! for Gultepe formula
     &  offset =(/     273d0,       595d0/)   ! for Gultepe formula

      integer it, n
      real*8  An,An0,cdnc(2),cdnc0(2),fbymass1

      do it=1,2  ! ocean, land
        An0 = anssdd(i,j)  !  aerosol number of sea salt and dust
        An  = An0          !  aerosol number of sea salt and dust
        do n=1,4
          fbymass1 =  F_act(n)*(.75d0/pi)/(dens(n)*radto3(it))
          An0 = An0 + mdpi (n,i,j)*fbymass1   ! +fact*tot_mass/part_mass
          An  = An  + mdcur(n,i,j)*fbymass1
        end do
        fbymass1 =  F_act(5)*(.75d0/pi)/(dens(5)*radto3(it))
        An  = An  + mdcur(5,i,j)*fbymass1

        if(An0.lt.1.) An0=1.
        if(An .lt.1.) An =1.
        cdnc0(it) = max( 20d0, scl(it)*log10(AN0)-offset(it))
        cdnc (it) = max( 20d0, scl(it)*log10(AN )-offset(it))
      end do

      dCDNC = (1-pland)*(cdnc(1)-cdnc0(1))+pland *(cdnc(2)-cdnc0(2))
      return
      end subroutine dCDNC_EST

      SUBROUTINE updateAerosol(JYEARA,JJDAYA, a6jday, plbaer)
      implicit none
      INTEGER, intent(in) :: jyeara,jjdaya
      real*8, pointer :: a6jday(:,:,:,:)
      real*8, dimension(:), pointer ::  plbaer

      call stop_model('updateAerosol: should not get here',255)

      RETURN
      END SUBROUTINE updateAerosol

      subroutine updateAerosol2(jYearA, jjDaya, a6jday, plbaer)
!@sum updateAerosol2 reads aerosol file(s) and calculates A6JDAY(lma,6,:,:)
!@+   (dry aerosol Tau) for current day, year.  On startup, it allocates
!@+   a6jday and plbaer, and reads plbaer.  Note that jYearA may be
!@+   negative, which is the Model E method for indicating that the
!@+   data for abs(jYearA) is to be used for all years (this matters
!@+   when time-interpolating between December and January).
! This version gives different results than the previous version, since
! the latter calculated the inter-month time interpolation weight using
!   XMI=(JJDAYA+JJDAYA+31-(JJDAYA+15)/61+(JJDAYA+14)/61)/61.D0
! and read_stream calculates it based on the midpoint dates of the
! current and following months.
! Note that the radiation code still assumes that all aerosols are on the
! same vertical grid.  A per-aerosol plbaer in calls to REPART in the
! radiation code would allow more flexibility.
      use domain_decomp_atm, only : grid
      use JulianCalendar_mod, only : jdmidofm ! for md1850 month interp
      use timestream_mod, only : init_stream,read_stream
     &     ,reset_stream_properties,get_by_index,getname_firstfile
      use pario, only : par_open,par_close,read_dist_data,read_data
     &     ,get_dimlen,get_dimlens
      implicit none
c
!@var jyeara current year; negative if the same year is to be repeated
!@var jjdaya current day
      INTEGER, intent(in) :: jyeara,jjdaya
!@var a6jday optical depth for 6 aerosol types
      real*8, pointer :: a6jday(:,:,:,:)
!@var plbaer pressures of layer interfaces in aerosol datafiles
      real*8, dimension(:), pointer ::  plbaer
c
      INTEGER m,mi,mj,i,j,l,n,jyearx
      REAL*8 wtmi,wtmj
      REAL*8 xsslt ! ,xdust
      real*8 :: dp,mindp
      logical, save :: init = .false.
      logical :: cyclic
      integer :: dlens(7)
      character(len=32) :: fname1_ssa

      integer :: i_0,i_1,j_0,j_1

      integer :: fid
      real*8, allocatable :: aerarr(:,:,:),arr12(:,:,:,:)

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      jyearx = abs(jyeara)
      
      if (.not. init) then
        init = .true.

        allocate( md1850 (4,i_0:i_1,j_0:j_1,0:12) )
        allocate(anssdd(i_0:i_1,j_0:j_1))
        allocate(mdpi(4,i_0:i_1,j_0:j_1))
        allocate(mdcur(5,i_0:i_1,j_0:j_1))

        cyclic = jyeara < 0

        ! Init the sea salt file first just to obtain lma,plbaer metadata
        n = 2
        call init_stream(grid,A6streams(n),
     &       'TAero_'//trim(aernames(n)),trim(aernames(n)),
     &       0d0,1d30,'linm2m',jyearx,jjdaya,cyclic=cyclic)
        call getname_firstfile(A6streams(n),fname1_ssa)

        fid = par_open(grid,trim(fname1_ssa),'read')

        !lma = get_dimlen(grid,fid,'lev')
        call get_dimlens(grid,fid,'plbaer',n,dlens)
        lma = dlens(1)-1

        if (.not.associated(A6JDAY))
     *       allocate(A6JDAY(lma,6,i_0:i_1,j_0:j_1))

        if (.not. associated(plbaer)) allocate( plbaer(lma+1) )
        call read_data(grid,fid,'plbaer',plbaer,bcast_all=.true.)

        call par_close(grid,fid)

!**** For parameterized AIE, find level whose top is closest to 720 mb
        mindp = 1d30
        do La720=1,lma
          dp = abs(plbaer(La720)-720d0)
          if(dp > mindp) exit
          mindp = dp
        enddo
        La720 = La720 - 2

        allocate(arr12(grid%i_strt_halo:grid%i_stop_halo,
     &                 grid%j_strt_halo:grid%j_stop_halo,lma,12))
        allocate(aerarr(i_0:i_1,j_0:j_1,12))

        do n=1,6
          if(n.eq.2) cycle ! skip sea salt
          ! Initialize the stream to the year 1850 to extract
          ! the monthly climatology of 1850 aerosols for parameterized AIE
          ! The read_stream call below will jump to the current year.
          call init_stream(grid,A6streams(n),
     &         'TAero_'//trim(aernames(n)),trim(aernames(n)),
     &         0d0,1d30,'linm2m',1850,1,cyclic=.true.)
          do m=1,12
            call get_by_index(grid,A6streams(n),m,arr12(:,:,:,m))
          enddo
          aerarr = byz_cm3 *
     &         SUM(arr12(i_0:i_1,j_0:j_1,1:La720,:), DIM=3)
          select case (n)
          case (1)
            md1850(1,:,:,1:12) = aerarr
          case (3,4,5)
            md1850(n-1,:,:,1:12) = aerarr
          case (6)
            md1850(4,:,:,1:12) = md1850(4,:,:,1:12) + aerarr
          end select
          ! Need this call to allow year jumps for cyclic case
          call reset_stream_properties(grid,A6streams(n),cyclic=cyclic)
        enddo

        md1850(:,:,:,0) = md1850(:,:,:,12)
        deallocate(arr12,aerarr)

      endif ! end init

C**** read and time-interpolate
      allocate(aerarr(grid%i_strt_halo:grid%i_stop_halo,
     &                grid%j_strt_halo:grid%j_stop_halo,lma))
      do n=1,6
        call read_stream(grid,A6streams(n),jyearx,jjdaya,aerarr)
        do j=j_0,j_1
        do i=i_0,i_1
        do l=1,lma
          a6jday(l,n,i,j)=aerarr(i,j,l)
        enddo
        enddo
        enddo
      enddo
      deallocate(aerarr)

! remove AERMIX scalings
      DO J=J_0,J_1
      DO I=I_0,I_1
      DO N=1,6
      DO L=1,lma
        A6JDAY(L,N,I,J)=(1000.D0*DRYM2G(N))*A6JDAY(L,N,I,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C**** Calculate terms for aerosol indirect effect parameterization

      do i=1,13
        if(jjdaya.le.jdmidofm(i)) then
          wtmi = real(jdmidofm(i)-jjdaya,kind=8)/
     &               (jdmidofm(i)-jdmidofm(i-1))
          wtmj = 1d0-wtmi
          mi = i-1
          if(mi > 11) mi=0
          mj = mi+1
          exit
        endif
      enddo

!!!   xdust=.33/(2000.*4.1888*(.40d-6)**3)     ! f/[rho*4pi/3*r^3] (/kg)
      xsslt=1.d0/(2000.*4.1888*(.44d-6)**3) ! x/particle-mass (/kg)

      do j=J_0,J_1
      do i=I_0,I_1

C**** sea salt contribution to anssdd
c SUM to L=5 for low clouds only
        anssdd(i,j) = SUM(A6JDAY(1:La720,2,I,J))/(1000.D0*DRYM2G(2))
     &       *byz_cm3 * Xsslt

C**** SU4,NO3,OCX,BCB,BCI (reordered: no sea salt, no pre-ind BCI)
        mdpi(:,i,j) =
     &       WTMI*md1850(:,i,j,mi) + WTMJ*md1850(:,i,j,mj)
        mdcur(1,i,j) = SUM (A6JDAY(1:La720,1,I,J))*
     &       byz_gcm3/drym2g(1)
        mdcur(2,i,j) = SUM (A6JDAY(1:La720,3,I,J))*
     &       byz_gcm3/drym2g(3)
        mdcur(3,i,j) = SUM (A6JDAY(1:La720,4,I,J))*
     &       byz_gcm3/drym2g(4)
        mdcur(4,i,j) = SUM (A6JDAY(1:La720,6,I,J))*
     &       byz_gcm3/drym2g(6)
        mdcur(5,i,j) = SUM (A6JDAY(1:La720,5,I,J))*
     &       byz_gcm3/drym2g(5)
      end do
      end do

C Misc. comments that were in the previous version of this routine
C                TROPOSPHERIC AEROSOL COMPOSITIONAL/TYPE PARAMETERS
C                   SO4    SEA    ANT    OCX    BCI    BCB   *BCB  *BCB
C     DATA REFDRY/0.200, 1.000, 0.300, 0.300, 0.100, 0.100, 0.200,0.050/
C
C     DATA REFWET/0.272, 1.808, 0.398, 0.318, 0.100, 0.100, 0.200,0.050/
C
C     DATA DRYM2G/4.667, 0.866, 4.448, 5.018, 9.000, 9.000, 5.521,8.169/
C
CKoch DATA DRYM2G/5.000, 2.866, 8.000, 8.000, 9.000, 9.000, 5.521,8.169/
C
C     DATA RHTMAG/1.788, 3.310, 1.756, 1.163, 1.000, 1.000, 1.000,1.000/
C
CRH70 DATA WETM2G/8.345, 2.866, 7.811, 5.836, 9.000, 9.000, 5.521,8.169/
C
C     DATA Q55DRY/2.191, 2.499, 3.069, 3.010, 1.560, 1.560, 1.914,0.708/
C
C     DATA DENAER/1.760, 2.165, 1.725, 1.500, 1.300, 1.300, 1.300,1.300/
C
C     ------------------------------------------------------------------
C          DRYM2G(I) = 0.75/DENAER(I)*Q55DRY(I)/REFDRY(I)
C          WETM2G(I) = DRYM2G(I)*RHTMAG(I)
C          RHTMAG(I) = Rel Humidity TAU Magnification factor  at RH=0.70
C          REFWET(I) = Rel Humidity REFDRY Magnification      at RH=0.70
C     ------------------------------------------------------------------

      return
      end subroutine updateAerosol2

      end module AerParam_mod

#ifdef OLD_BCdalbsn
      subroutine updBCd(year)
!@sum updBCd reads timeseries file for black carbon deposition
!@+   and interpolates depoBC to requested year.
!@auth R. Ruedy, M. Kelley
      use domain_decomp_atm, only : grid,getDomainBounds
      use timestream_mod, only : init_stream,read_stream
      use AerParam_mod, only: BCdepstream,depoBC,depoBC_1990
      implicit none
      integer, intent(in) :: year
c
      logical, save :: init = .false.
      integer :: i_0h,i_1h,j_0h,j_1h
      integer :: day

      day = 1 ! to pass a required argument

      if (.not. init) then
        init = .true.

        call getDomainBounds(grid, i_strt_halo=i_0h, i_stop_halo=i_1h,
     &                             j_strt_halo=j_0h, j_stop_halo=j_1h)
        allocate(depoBC     (i_0h:i_1h, j_0h:j_1h),
     &           depoBC_1990(i_0h:i_1h, j_0h:j_1h))

        call init_stream(grid,BCdepstream,'BC_dep','BC_dep',
     &       0d0,1d30,'none',year,day)
      endif

      call read_stream(grid,BCdepstream,year,day,depoBC)

      end subroutine updBCd
#else
      subroutine updBCdalbsn(year,day)
!@sum updBCdalbsn reads timeseries file for black carbon delta-snow-albedo
!@+   and interpolates to the requested day/year.   If the year is negative,
!@+   this is interpreted as indicating perpetual-year mode, as per the
!@+   convention for numerous radiation input files.
!@auth R. Ruedy, M. Kelley
      use domain_decomp_atm, only : grid,getDomainBounds
      use timestream_mod, only : init_stream,read_stream
      use AerParam_mod, only: BCdalbsnstream,BCdalbsn
      implicit none
      integer, intent(in) :: year,day
c
      logical, save :: init = .false.
      integer :: i_0h,i_1h,j_0h,j_1h
      logical :: cyclic
      integer :: absyr

      absyr = abs(year)
      if (.not. init) then
        init = .true.

        call getDomainBounds(grid, i_strt_halo=i_0h, i_stop_halo=i_1h,
     &                             j_strt_halo=j_0h, j_stop_halo=j_1h)
        allocate(BCdalbsn(i_0h:i_1h, j_0h:j_1h))
        BCdalbsn = 0.
        cyclic = year < 0
        call init_stream(grid,BCdalbsnstream,'BCdalbsn','BCdalbsn',
     &       -1d30,1d30,'linm2m',absyr,day,cyclic=cyclic)
      endif

      call read_stream(grid,BCdalbsnstream,absyr,day,BCdalbsn)
      BCdalbsn = BCdalbsn / 100d0 ! units conversion from % to 1

      end subroutine updBCdalbsn
#endif

      module DustParam_mod
!@sum This module reads, time-interpolates, and stores fields needed
!@+   by the radiation code in the prescribed-dust configuration
!@+   of modelE.   The logic follows that of AerParam_mod.
!@+   The interface routine is upddst2().
!@auth R. Miller original version
!@auth M. Kelley reprogrammed for new-style time-varying input
      use timestream_mod, only : timestream
      implicit none

!@var ddjday (kg/m2/layer) dust amount for each size class, layer, and column
!@+   for the current day
      real*8, dimension(:,:,:,:), allocatable :: ddjday

!@var {lmd,nsized} number of {layers, size classes} in DUSTaer input file
      integer :: lmd,nsized

!@var DUSTaerstream interface for reading and time-interpolating DUSTaer files
!@+   See usage notes in timestream_mod
      type(timestream) :: DUSTaerstream

!@var is_initialized whether the DUSTaer stream has been initialized
!@+   and various arrays allocated
      logical :: is_initialized=.false.

!@var plbdust nominal edge pressures of DUSTaer file layers
!@var {re,ro}dust radii of DUSTaer file size classes
      real*8, dimension(:), allocatable :: redust, rodust, plbdust

      contains

      subroutine upddst2(jyeard,jjdayd)
      use domain_decomp_atm, only : grid
      use timestream_mod, only : init_stream,read_stream,
     &     getname_firstfile
      use pario, only : par_open,par_close
     &     ,get_dimlens,read_data
      implicit none
!@var jyeard, jjdayd year and day of the data to read into ddjday.
!@+   Note that jyeard may be negative, which is the Model E method for
!@+   indicating that the data for abs(jyeard) is to be used for all years.
      integer, intent(in) :: jyeard,jjdayd
!
      integer :: i_0,i_1,j_0,j_1,i,j,n,jyearx
      logical :: cyclic
      real*8, dimension(:,:,:,:), allocatable :: ddjday_transp
      integer :: fid,ndims,dlens(7)
      character(len=32) :: fname1_dust

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      cyclic = jyeard < 0
      jyearx = abs(jyeard)

      if(.not. is_initialized) then
        is_initialized = .true.

        call init_stream(grid,DUSTaerstream,
     &       'DUSTaer','DUST',
     &       0d0,1d30,'linm2m',jyearx,jjdayd,cyclic=cyclic)


        ! read dust metadata
        call getname_firstfile(DUSTaerstream,fname1_dust)
        fid = par_open(grid,trim(fname1_dust),'read')
        call get_dimlens(grid,fid,'DUST',ndims,dlens)
        lmd    = dlens(3)
        nsized = dlens(4)
        allocate( plbdust(lmd+1), redust(nsized), rodust(nsized) )
        call read_data(grid,fid,'plbdust',plbdust,bcast_all=.true.)
        call read_data(grid,fid,'redust',redust,bcast_all=.true.)
        call read_data(grid,fid,'rodust',rodust,bcast_all=.true.)
        call par_close(grid,fid)

        allocate( ddjday(lmd,nsized,i_0:i_1,j_0:j_1) )

      endif

      allocate( ddjday_transp(
     &     grid%i_strt_halo:grid%i_stop_halo,
     &     grid%j_strt_halo:grid%j_stop_halo,
     &     lmd,nsized) )

      call read_stream(grid,DUSTaerstream,jyearx,jjdayd,ddjday_transp)

      do j=j_0,j_1
      do i=i_0,i_1
        ddjday(:,:,i,j) = ddjday_transp(i,j,:,:)
      enddo
      enddo

      deallocate ( ddjday_transp )

      end subroutine upddst2

      end module DustParam_mod
