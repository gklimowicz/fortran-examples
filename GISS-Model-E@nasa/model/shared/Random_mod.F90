module RANDOM
!@sum   RANDOM generates random numbers: 0<RANDom_nUmber<1
!@auth  Reto Ruedy
!@ver   1.0 (SGI,IBM,Linux,DEC)
!@cont  RANDU, RINIT, RFINAL
  implicit none
  integer, save :: IX            !@var IX     random number seed

  ! Parameters used for generating sequences of random numbers
  integer, parameter :: LONG = selected_int_kind(12)
  integer(kind=LONG), parameter :: A_linear = 69069 ! Based on ifort ran() routine
  integer, parameter :: MAX_BITS = 32
  real*4, parameter :: SCALE = 1./real(2_8**MAX_BITS,4)

contains

  real*8 function RANDU (x)
!@sum   RANDU calculates a random number based on the seed IX
    real*8 :: x ! unused
    integer, parameter :: MASK = int(z'ffffff00')

    ix = ix * A_LINEAR + 1

    if (ix < 0) then
      RANDU = 1.0 + (iand(ix, MASK) * SCALE)
    else
      RANDU = iand(ix, MASK) * SCALE
    end if

  end function RANDU

  subroutine RINIT (INIT)
!@sum   RINIT sets the initial seed IX
    integer, intent(IN)  :: INIT   !@var INIT   first random no. seed
    IX=INIT
    return
  end subroutine RINIT

  subroutine RFINAL (IFINAL)
!@sum   RFINAL retrieves final seed value
    integer, intent(OUT) :: IFINAL !@var IFINAL last random no. seed
    IFINAL=IX
    return
  end subroutine RFINAL

  subroutine BURN_RANDOM(n)
!@sum  BURN_RANDOM burns a set number of random numbers. It is used to
    !                  maintain bit-wise correspondence on parallel runs.
    implicit none
    integer, intent(in) :: n
    integer :: i
    real*8 x, randss
    integer :: a, b ! linear coefficient
    integer :: nn
    if (n.eq.0) return
    a = A_linear
    nn = n
    b = 1
    do i=1,MAX_BITS
      if (mod(nn,2) == 1) ix = ix * a + b
      b=(a+1)*b
      a=a*a
      nn=nn/2
      if (nn == 0) exit
    end do
    return
  end subroutine burn_random

end module RANDOM
