module Precision_mod
!@sum  The reduce_precision routines truncate the number of
!@+    significant digits in a real*8 number x to an approximate
!@+    precision of relacc (1d-16 < relacc << 1).  Fortran functions
!@+    can be used to define x = fraction(x) * 2**exponent(x),
!@+    where fraction(x) = O(1).  The part of fraction(x) smaller than
!@+    relacc is discarded.
!@auth M. Kelley
  implicit none
  public :: reduce_precision
  interface reduce_precision
    module procedure reduce_precision_0d
    module procedure reduce_precision_1d
    module procedure reduce_precision_2d
    module procedure reduce_precision_3d
    module procedure reduce_precision_4d
  end interface
contains
  subroutine reduce_precision_0d(x,relacc)
    real*8, intent(inout) :: x
    real*8, intent(in) :: relacc
    x = nint(fraction(x)/relacc,kind=8)*relacc*2d0**exponent(x)
  end subroutine reduce_precision_0d
  subroutine reduce_precision_1d(x,relacc)
    real*8, intent(inout) :: x(:)
    real*8, intent(in) :: relacc
    x = nint(fraction(x)/relacc,kind=8)*relacc*2d0**exponent(x)
  end subroutine reduce_precision_1d
  subroutine reduce_precision_2d(x,relacc)
    real*8, intent(inout) :: x(:,:)
    real*8, intent(in) :: relacc
    x = nint(fraction(x)/relacc,kind=8)*relacc*2d0**exponent(x)
  end subroutine reduce_precision_2d
  subroutine reduce_precision_3d(x,relacc)
    real*8, intent(inout) :: x(:,:,:)
    real*8, intent(in) :: relacc
    x = nint(fraction(x)/relacc,kind=8)*relacc*2d0**exponent(x)
  end subroutine reduce_precision_3d
  subroutine reduce_precision_4d(x,relacc)
    real*8, intent(inout) :: x(:,:,:,:)
    real*8, intent(in) :: relacc
    x = nint(fraction(x)/relacc,kind=8)*relacc*2d0**exponent(x)
  end subroutine reduce_precision_4d

end module Precision_mod
