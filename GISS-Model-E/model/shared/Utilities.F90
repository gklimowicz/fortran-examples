#include "rundeck_opts.h"
!@sum  Model Independent Utilities
!@auth Original Development Team

function THBAR (X,Y)
!@sum  THBAR calculates mean temperature used in vertical differencing
!@auth Gary Russell, Jean Lerner, Arakawa
  !****
  !**** THBAR(T1,T2) = (ln(T1) - ln(T2))/(1/T2 - 1/T1)
  !****              = T1*g(x) with x=T1/T2 , g(x)=ln(x)/(x-1)
  !****      g(x) is replaced by a rational function
  !****           (a+bx+cxx+dxxx+xxxx)/(e+fx+gxx)
  !****      approx.error <1.E-6 for x between .9 and 1.7
  !****
  implicit none
!@var A,B,C,D,E,F,G   expansion coefficients for THBAR
  real*8, parameter :: A=113.4977618974100d0
  real*8, parameter :: B=438.5012518098521d0
  real*8, parameter :: C=88.49964112645850d0
  real*8, parameter :: D=-11.50111432385882d0
  real*8, parameter :: E=30.00033943846368d0
  real*8, parameter :: F=299.9975118132485d0
  real*8, parameter :: G=299.9994728900967d0
  real*8 :: Q,AL                 !@var Q,AL   working variables
  real*8, intent(IN) :: X,Y      !@var X,Y    input temperatures
  real*8 :: THBAR                !@var THBAR  averaged temperature
  Q=X/Y
  AL=(A+Q*(B+Q*(C+Q*(D+Q))))/(E+Q*(F+G*Q))
  THBAR=X*AL
  return
end function THBAR

function QSAT (TM,LH,PR)
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
  use CONSTANT, only : mrat,rvap,tf
  implicit none
!@var A,B,C   expansion coefficients for QSAT
  real*8, parameter :: A=6.108d0*MRAT    !3.797915d0
  real*8, parameter :: B= 1./(RVAP*TF)   !7.93252d-6
  real*8, parameter :: C= 1./RVAP        !2.166847d-3
  !**** Note that if LH is considered to be a function of temperature, the
  !**** correct argument in QSAT is the average LH from t=0 (C) to TM, ie.
  !**** LH = 0.5*(LH(0)+LH(t))
  real*8, intent(IN) :: TM  !@var TM   temperature (K)
  real*8, intent(IN) :: PR  !@var PR   air pressure (mb)
  real*8, intent(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
  real*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
  QSAT = A*exp(LH*(B-C/max(130.d0,TM)))/PR
  return
end function QSAT

function DQSATDT (TM,LH)
!@sum  DQSATDT calculates change of sat. vapour mixing ratio with temp.
!@auth Gary Russell
  !**** Note that d(qsat)/dt = qsat * lh * c / T*T
  !**** Only the factor of qsat is given here
  use CONSTANT, only : rvap
  implicit none
!@var C coefficient for QSAT
  real*8, parameter :: C = 1./RVAP        !2.166847d-3
  !**** Note that if LH is considered to be a function of temperature, the
  !**** correct argument in DQSATDT is the actual LH at TM i.e. LH=LH(TM)
  real*8, intent(IN) :: TM  !@var TM   temperature (K)
  real*8, intent(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
  real*8 :: DQSATDT         !@var DQSATDT d(qsat)/dT factor only.
  DQSATDT = LH*C/(TM*TM)    ! * QSAT(TM,LH,PR)
  return
end function DQSATDT

function wv_psat( tm, lh )
!@sum wv_psat calculates saturation water vapor pressure
!@auth Gary Russell (qsat), Jan Perlwitz

  use constant, only : rvap,tf

  implicit none

!@var A,B,C   expansion coefficients
  real( kind=8 ), parameter :: A = 6.108d0      ![hPa]
  real( kind=8 ), parameter :: B = 1./(RVAP*TF) !7.93252d-6
  real( kind=8 ), parameter :: C = 1./RVAP      !2.166847d-3
!**** Note that if LH is considered to be a function of temperature, the
!**** correct argument in wv_psat is the average LH from t=0 (C) to TM, ie.
!**** LH = 0.5*(LH(0)+LH(t))
!@var TM   temperature (K)
  real( kind=8 ), intent(IN) :: TM
!@var LH   lat. heat of vap./sub. (J/kg)
  real( kind=8 ), intent(IN) :: LH
!@var wv_psat saturation water vapor pressure [hPa]
  real( kind=8 ) :: wv_psat

  wv_psat = A*exp(LH*(B-C/max(130.d0,TM)))

  return
end function wv_psat

function SLP(PS,TAS,ZS)
!@sum SLP estimates sea level pressure in the presence of topography
!@+   for better match to reanalyses.
  use CONSTANT, only: bmoist, grav, rgas, by3
  implicit none
!@var PS surface pressure (mb)
!@var TAS surface temperature (K)
!@var ZS surface elevation (m)
  real*8, intent(IN) :: PS, TAS, ZS
  real*8 :: SLP, TSL, BETA, BZBYT, GBYRB, TASn

  if (ZS.ne.0.) then
    TSL= TAS+BMOIST*ZS
    TASn=TAS
    BETA=BMOIST
    if (TAS < 290.5 .and. TSL > 290.5) BETA= (290.5d0 - TAS)/ZS
    if (TAS > 290.5 .and. TSL > 290.5) TASn = 0.5*(290.5d0 + TAS)
    if (TAS < 255) TASn = 0.5*(255d0 + TAS)
    BZBYT=BETA*ZS/TASn
    GBYRB=GRAV/(RGAS*BETA)
    if (BETA > 1d-6 ) then
      SLP=PS*(1.+BZBYT)**GBYRB
    else
      SLP=PS*exp((1.-0.5*BZBYT+BZBYT**(2.*by3))*GBYRB*BZBYT)
    end if
  else
    SLP=PS
  end if
  return
end function SLP


subroutine READT (IUNIT,NSKIP,LENGTH,AOUT,IPOS)
!@sum   READT  read in title and real*4 array and convert to real*8
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
  use FILEMANAGER, only : NAME=>nameunit
  implicit none
  integer, intent(IN) :: IUNIT        !@var  IUNIT  file unit number
  integer, intent(IN) :: NSKIP    !@var  NSKIP  no. of R*4's to skip
  integer, intent(IN) :: LENGTH       !@var  LENGTH size of array
  integer, intent(IN) :: IPOS  !@var  IPOS   no. of recs. to advance
  !REAL*4, INTENT(OUT) :: AIN(LENGTH)  !@var  AIN    real*4 array
  real*8, intent(OUT) :: AOUT(LENGTH) !@var  AOUT   real*8 array
  real*4 :: X               !@var  X      dummy variable
  integer :: N              !@var  N      loop variable
  character*80 TITLE        !@var  TITLE  title of file record
  real*4, allocatable :: buf(:)

  allocate( buf(LENGTH) )

  do N=1,IPOS-1
    read (IUNIT,end=920)
  end do
  read (IUNIT,ERR=910,end=920) TITLE,(X,N=1,NSKIP),buf
  !**** do transfer backwards in case AOUT and AIN are same workspace
!!! THIS DOESN''T WORK IN F90+  !!! stop using such hacks !
  !      DO N=LENGTH,1,-1
  !        AOUT(N)=AIN(N)
  !      END DO
  AOUT(:) = buf(:)
  deallocate( buf )
  write(6,*) "Read from file ",trim(NAME(IUNIT)),": ",trim(TITLE)
  return
910 write(6,*) 'READ ERROR ON FILE ',NAME(IUNIT)
  call stop_model('tREAD: READ ERROR',255)
920 write(6,*) 'END OF FILE ENCOUNTERED ON FILE ',NAME(IUNIT)
  call stop_model('tREAD: No data found',255)
end subroutine READT

subroutine WRITEI (iunit,it,aout,len4)
!@sum   WRITEI  writes array surrounded by IT and secures it
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
  use FILEMANAGER, only : NAME=>nameunit
  implicit none
  integer, intent(IN) :: IUNIT       !@var  IUNIT  file unit number
  integer, intent(IN) :: IT          !@var  IT time, 1st & last word
  integer, intent(IN) :: LEN4        !@var  LENGTH size of array
  real*4,  intent(IN) :: AOUT(LEN4)  !@var  AOUT   real*4 array

  write (iunit) it,aout,it
  call sys_flush(iunit)
  write (6,*) "Wrote to file ",trim(NAME(IUNIT)),", time=",it
  return
end subroutine WRITEI

subroutine READI (iunit,it,ain,it1,len4,iok)
!@sum  READI reads array surrounded by IT (for post processing)
!@auth  Original Development Team
!@ver   1.0
  implicit none
  integer, intent(IN) :: IUNIT,LEN4
  integer, intent(out) :: IT,it1,iok ! iok: 0='ok',1='not ok'
  real*4  AIN(LEN4)
  iok = 0
  read(iunit,end=555) it,ain,it1
  return
555 iok=1
  return
end subroutine readi

subroutine WRITEI8 (iunit,it,aout,len8)
!@sum   WRITEI8 writes real*8 array surrounded by IT and secures it
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
  use FILEMANAGER, only : NAME=>nameunit
  implicit none
  integer, intent(IN) :: IUNIT       !@var  IUNIT  file unit number
  integer, intent(IN) :: IT          !@var  IT time, 1st & last word
  integer, intent(IN) :: LEN8        !@var  LENGTH size of array
  real*8,  intent(IN) :: AOUT(LEN8)  !@var  AOUT   real*8 array

  write (iunit) it,aout,it
  call sys_flush(iunit)
  write (6,*) "Wrote to file ",trim(NAME(IUNIT)),", time=",it
  return
end subroutine WRITEI8

subroutine READI8 (iunit,it,ain,it1,len8,iok)
!@sum  READI reads array surrounded by IT (for post processing)
!@auth  Original Development Team
!@ver   1.0
  implicit none
  integer, intent(IN) :: IUNIT,LEN8
  integer, intent(out) :: IT,it1,iok ! iok: 0='ok',1='not ok'
  real*8  AIN(LEN8)
  iok = 0
  read(iunit,end=555) it,ain,it1
  return
555 iok=1
  return
end subroutine readi8

subroutine io_POS (iunit,it,len4,itdif)
!@sum   io_POS  positions a seq. output file for the next write operat'n
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
  use FILEMANAGER, only : NAME=>nameunit
  implicit none
  integer, intent(IN) :: IUNIT !@var IUNIT  file unit number
  integer, intent(IN) :: IT,ITdif !@var IT,ITdif current time,dt
  integer, intent(IN) :: LEN4 !@var LENGTH of array in words
  integer :: IT1,IT2   !@var  time_tags at start,end of each record
  integer :: N              !@var  N      loop variable

  read (iunit,end=10,err=50) it1,(it2,n=1,len4+1)
  if(it1 .le. it) go to 30
10 write(6,*) "Starting a new file ",trim(NAME(IUNIT)),", time=",it
  rewind iunit
  return

20 read (iunit,end=35,err=50) it1,(it2,n=1,len4+1)
30 if (it2 .ne. it1) then
    write(6,*) 'file ',trim(NAME(IUNIT)),' damaged: it/it1/it2=', &
         &    it,it1,it2
    call stop_model('io_POS: damaged file',255)
  end if
  if (it1 .le. it) go to 20
  it1=it1-itdif
35 backspace iunit
  if (it1+itdif .le. it) go to 40
  write (6,*) "positioned ",trim(NAME(IUNIT)),", it1/itime=",it1,it
  return
40 write (6,*) "file ",trim(NAME(IUNIT))," too short, it1/it=",it1,it
  call stop_model('io_POS: file too short',255)
50 write (6,*) "Read error on: ",trim(NAME(IUNIT)),", it1/it=",it1,it
  call stop_model('io_POS: read error',255)
end subroutine io_POS


subroutine CHECK3(A,IN,JN,LN,SUBR,FIELD)
!@sum  CHECK3 Checks for NaN/INF in real 3-D arrays
!@auth Original development team
  implicit none

!@var IN,JN,LN size of 3-D array
  integer, intent(IN) :: IN,JN,LN
!@var SUBR identifies where CHECK3 was called from
  character*(*), intent(IN) :: SUBR
!@var FIELD identifies the field being tested
  character*(*), intent(IN) :: FIELD
!@var A array being tested
  real*8, dimension(IN,JN,LN),intent(IN) :: A
  logical :: QCHECK3 = .false.
  integer I,J,L !@var I,J,L loop variables

  do L=1,LN
    do J=1,JN
      do I=1,IN
        if (.not.(A(I,J,L).gt.0..or.A(I,J,L).le.0.) .or. &
             &       abs(A(I,J,L)) .gt.huge(A(I,J,L)) ) then
          Write (6,9) Trim(FIELD),I,J,L,A(I,J,L),Trim(SUBR)
          if (J.lt.JN.and.J.gt.1) QCHECK3 = .true.
        end if
      end do
    end do
  end do
  call SYS_FLUSH(6)
  if (QCHECK3) call stop_model('CHECK3',255)
  return
9 Format ('CHECK3: FIELD,I,J,L,VALUE = ',A,3I6,E16.8,'  after  ',A)
end subroutine CHECK3


subroutine CHECK3B(A,I1,I2,J1,J2,NJPOL,LN,SUBR,FIELD)
!@sum  CHECK3B Checks for NaN/INF in real 3-D arrays
!@auth Original development team
  implicit none

!@var IN,JN,LN size of 3-D array
  integer, intent(IN) :: I1,I2,J1,J2,NJPOL,LN
!@var SUBR identifies where CHECK3 was called from
  character*(*), intent(IN) :: SUBR
!@var FIELD identifies the field being tested
  character*(*), intent(IN) :: FIELD
!@var A array being tested
  real*8, dimension(I1:I2,J1:J2,LN),intent(IN) :: A
  logical :: QCHECK3 = .false.
  integer I,J,L !@var I,J,L loop variables

  do L=1,LN
    do J=J1+NJPOL,J2-NJPOL
      do I=I1,I2
        if (.not.(A(I,J,L).gt.0..or.A(I,J,L).le.0.) .or. &
             &       abs(A(I,J,L)) .gt.huge(A(I,J,L)) ) then
          Write (6,9) Trim(FIELD),I,J,L,A(I,J,L),Trim(SUBR)
          QCHECK3 = .true.
        end if
      end do
    end do
  end do
  call SYS_FLUSH(6)
  if (QCHECK3) call stop_model('CHECK3B',255)
  return
9 Format ('CHECK3B: FIELD,I,J,L,VALUE = ',A,3I6,E16.8,'  after  ',A)
end subroutine CHECK3B


subroutine CHECK3C(A,LN,I1,I2,J1,J2,NJPOL,SUBR,FIELD)
!@sum  CHECK3B Checks for NaN/INF in real 3-D arrays
!@auth Original development team
  implicit none

!@var IN,JN,LN size of 3-D array
  integer, intent(IN) :: LN,I1,I2,J1,J2,NJPOL
!@var SUBR identifies where CHECK3 was called from
  character*(*), intent(IN) :: SUBR
!@var FIELD identifies the field being tested
  character*(*), intent(IN) :: FIELD
!@var A array being tested
  real*8, dimension(LN,I1:I2,J1:J2),intent(IN) :: A
  logical :: QCHECK3 = .false.
  integer I,J,L !@var I,J,L loop variables

  do J=J1+NJPOL,J2-NJPOL
    do I=I1,I2
      do L=1,LN
        if (.not.(A(L,I,J).gt.0..or.A(L,I,J).le.0.) .or. &
             &       abs(A(L,I,J)) .gt.huge(A(L,I,J)) ) then
          Write (6,9) Trim(FIELD),L,I,J,A(L,I,J),Trim(SUBR)
          QCHECK3 = .true.
        end if
      end do
    end do
  end do
  call SYS_FLUSH(6)
  if (QCHECK3) call stop_model('CHECK3C',255)
  return
9 Format ('CHECK3C: FIELD,L,I,J,VALUE = ',A,3I6,E16.8,'  after  ',A)
end subroutine CHECK3C


subroutine CHECK4(A,IN,JN,KN,LN,SUBR,FIELD)
!@sum  CHECK4 Checks for NaN/INF in real 4-D arrays
!@auth Original development team
  implicit none

!@var IN,JN,KN,LN size of 4-D array
  integer, intent(IN) :: IN,JN,KN,LN
!@var SUBR identifies where CHECK4 was called from
  character*(*), intent(IN) :: SUBR
!@var FIELD identifies the field being tested
  character*(*), intent(IN) :: FIELD
!@var A array being tested
  real*8, dimension(IN,JN,KN,LN),intent(IN) :: A
  logical :: QCHECK4 = .false.
  integer I,J,K,L !@var I,J,K,L loop variables

  do L=1,LN
    do K=1,KN
      do J=1,JN
        do I=1,IN
          if (.not.(A(I,J,K,L).gt.0..or.A(I,J,K,L).le.0.) .or. &
               &       abs(A(I,J,K,L)) .gt.huge(A(I,J,K,L)) ) then
            Write (6,*) Trim(FIELD),I,J,K,L,A(I,J,K,L),Trim(SUBR)
            if (J.lt.JN.and.J.gt.1) QCHECK4 = .true.
          end if
        end do
      end do
    end do
  end do
  call SYS_FLUSH(6)
  if (QCHECK4) call stop_model('CHECK4',255)
  return
9 Format ('CHECK4: FIELD,I,J,K,L,VALUE = ',A,3I6,E16.8,'  after  ',A)
end subroutine CHECK4


subroutine CHECK4B(A,I1,I2,J1,J2,NJPOL,KN,LN,SUBR,FIELD)
!@sum  CHECK4 Checks for NaN/INF in real 4-D arrays
!@auth Original development team
  implicit none

!@var IN,JN,KN,LN size of 4-D array
  integer, intent(IN) :: I1,I2,J1,J2,NJPOL,KN,LN
!@var SUBR identifies where CHECK4 was called from
  character*(*), intent(IN) :: SUBR
!@var FIELD identifies the field being tested
  character*(*), intent(IN) :: FIELD
!@var A array being tested
  real*8, dimension(I1:I2,J1:J2,KN,LN),intent(IN) :: A
  logical :: QCHECK4 = .false.
  integer I,J,K,L !@var I,J,K,L loop variables

  do L=1,LN
    do K=1,KN
      do J=J1+NJPOL,J2-NJPOL
        do I=I1,I2
          if (.not.(A(I,J,K,L).gt.0..or.A(I,J,K,L).le.0.) .or. &
               &       abs(A(I,J,K,L)) .gt.huge(A(I,J,K,L)) ) then
            Write (6,*) Trim(FIELD),I,J,K,L,A(I,J,K,L),Trim(SUBR)
            QCHECK4 = .true.
          end if
        end do
      end do
    end do
  end do
  call SYS_FLUSH(6)
  if (QCHECK4) call stop_model('CHECK4B',255)
  return
9 Format ('CHECK4B: FIELD,I,J,L,VALUE = ',A,3I6,E16.8,'  after  ',A)
end subroutine CHECK4B


function unit_string (pow10,ending)
!@sum Construct a units string with nice properties (no embedded blanks)
!@auth G. Schmidt, J. Lerner
  !**** If a trailing ')' is supplied, it is assumed that a leading
  !****      '(' is required, so it is inserted
  implicit none
  character*(*) ending,unit_string
  character*10 tpow
  integer pow10,p

  tpow = ' '
  if(pow10.ne.0) then
    write(tpow,'(i3)') pow10
    tpow= '10^'//trim(adjustl(tpow))
    p=len_trim(ending)
    if (p > 0) then
      if (ending(p:p)==')') then
        tpow='('//trim(adjustl(tpow))
      end if
    end if
  endif
  unit_string = adjustl(trim(tpow)//" "//trim(ending))
  return
end function unit_string



subroutine write_run_status( message, retcode )
  implicit none
  character*(*), intent (in) :: message
  integer, intent(in) :: retcode
  integer, parameter :: iu_status = 9
  character*10 :: form_str
  integer num_digits

  ! construct format string in such a way that retcode is printed
  ! at the beginning of the line with no extra spaces
  if ( retcode .ne. 0 ) then
    num_digits = log10( real( abs(retcode), kind(1.d0) ) ) + 1
  else
    num_digits = 1
  endif
  if ( retcode < 0 ) num_digits = num_digits + 1

  write(form_str,"('(I',I1,')')") num_digits

  open( iu_status, file='run_status', form='FORMATTED', &
       &     status='UNKNOWN', ERR=10 )
  write( iu_status, form_str, ERR=10 ) retcode
  write( iu_status, '(A)', ERR=10 ) message
  close( iu_status )

  return
10 continue
  write( 0, * ) "ERROR: Can't write to the run_status file"
  write( 0, * ) "STATUS:", message
end subroutine write_run_status


function clean_str(string)
!@sum clean_str utility to clean strings of netcdf-unfriendly characters
  implicit none
  character(len=*), intent(in) :: string
  character(len=40) :: clean_str
  integer :: k

  clean_str=trim(string)
  do k=1,len_trim(string)
    if (clean_str(k:k).eq." " .or. clean_str(k:k).eq."+") &
         &        clean_str(k:k)="_"
  end do
  return
end function clean_str

function i5toc4 (iyr)
!@sum i5toc4 converts an integer 0-35999 into a character string of length 4
!@+   to preserve output file naming convention if runs go past year 9999
!@auth Reto Ruedy
      integer, intent(in) :: iyr
      character(len=4) i5toc4
      integer iyrbyk
      character(len=26) :: atoz='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

      if(iyr<10000) then
         write(i5toc4,'(i4.4)') iyr
      else if(iyr>35999) then
         i5toc4='xxxx'
      else
         iyrbyk=iyr/1000 - 9 ; i5toc4(1:1)=atoz(iyrbyk:iyrbyk)
         write(i5toc4(2:4),'(i3.3)') mod(iyr,1000)
      end if

      return
end function i5toc4


