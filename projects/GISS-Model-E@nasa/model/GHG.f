      module ghg_mod
      implicit none
      private

      public GTREND, ghg_init, CO2_trend

      INTEGER, PARAMETER :: nghg=6
      INTEGER ghgyr1,ghgyr2
      real*8, allocatable :: ghgam(:,:)

      contains

      subroutine CO2_trend(xCO2, TNOW)
      REAL*8, intent(out)  :: xCO2
      REAL*8, intent(in)  :: TNOW
      real*8 :: XNOW(nghg)
      
      call GTREND(XNOW,TNOW)
      xCO2 = XNOW(1)

      end subroutine CO2_trend


      SUBROUTINE GTREND(XNOW,TNOW)
C
!!      USE RADPAR, only: nghg,ghgyr1,ghgyr2,ghgam
      REAL*8 xnow(nghg),tnow,year,dy,frac
      INTEGER iy,n
      logical, save :: initialized=.false.

      if ( .not. initialized ) then
        call ghg_init
        initialized = .true.
      endif
C
C-------------------------------------------------------------
C        Makiko GHG Trend Compilation  GHG.1850-2050.Dec1999
C
C        Annual-Mean      Greenhouse Gas Mixing Ratios
C-------------------------------------------------------------
C                 CO2     N2O     CH4   CFC-11  CFC-12  others
C        Year     ppm     ppm     ppm     ppb     ppb     ppb
C-------------------------------------------------------------
C     Read from external file - outside table: use value from
C                                      years ghgyr1 or ghgyr2
      YEAR=TNOW
      IF(TNOW <= ghgyr1+.5D0) YEAR=ghgyr1+.5D0
      IF(TNOW >= ghgyr2+.49999D0) YEAR=ghgyr2+.49999D0
      DY=YEAR-(ghgyr1+.5D0)
      IY=DY
      frac=DY-IY
      IY=IY+1
C
C     CO2 N2O CH4 CFC-11 CFC-12 other_GHG  SCENARIO
C--------------------------------------------------
C
      do n=1,nghg
        XNOW(N)=GHGAM(N,IY)+frac*(GHGAM(N,IY+1)-GHGAM(N,IY))
      end do
C
      RETURN
      END SUBROUTINE GTREND


      SUBROUTINE GHGHST(iu)
!@sum  reads history for nghg well-mixed greenhouse gases
!@auth R. Ruedy

      use domain_decomp_atm, only : write_parallel
!!      USE RADPAR, only : nghg,ghgyr1,ghgyr2,ghgam
cddd      USE RAD_COM, only : ghg_yr
      INTEGER :: iu,n,k,nhead=4,iyr
      CHARACTER*80 title
      character(len=300) :: out_line

      write(out_line,*)  ! print header lines and first data line
      call write_parallel(trim(out_line),unit=6)
      do n=1,nhead+1
        read(iu,'(a)') title
        write(out_line,'(1x,a80)') title
        call write_parallel(trim(out_line),unit=6)
      end do
      if(title(1:2).eq.'--') then                 ! older format
        read(iu,'(a)') title
        write(out_line,'(1x,a80)') title
        call write_parallel(trim(out_line),unit=6)
        nhead=5
      end if

!**** find range of table: ghgyr1 - ghgyr2
      read(title,*) ghgyr1
      do ; read(iu,'(a)',end=20) title ; end do
   20 read(title,*) ghgyr2
      rewind iu  !   position to data lines
      do n=1,nhead ; read(iu,'(a)') ; end do

      allocate (ghgam(nghg,ghgyr2-ghgyr1+1))
      do n=1,ghgyr2-ghgyr1+1
        read(iu,*) iyr,(ghgam(k,n),k=1,nghg)
        do k=1,nghg ! replace -999. by reasonable numbers
          if(ghgam(k,n).lt.0.) ghgam(k,n)=ghgam(k,n-1)
        end do
cddd        if(ghg_yr>0 .and. abs(ghg_yr-iyr).le.1) then
cddd          write(out_line,'(i5,6f10.4)') iyr,(ghgam(k,n),k=1,nghg)
cddd          call write_parallel(trim(out_line),unit=6)
cddd        endif
      end do
      write(out_line,*) 'read GHG table for years',ghgyr1,' - ',ghgyr2
      call write_parallel(trim(out_line),unit=6)
      return
      end SUBROUTINE GHGHST


      subroutine ghg_init
      USE FILEMANAGER, only: openunit,closeunit
      integer :: iu

      call openunit('GHG',iu,.false.,.true.)
      call ghghst(iu)
      call closeunit(iu)

      end subroutine ghg_init

      end module ghg_mod
