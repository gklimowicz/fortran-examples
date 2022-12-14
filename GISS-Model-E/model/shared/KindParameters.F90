module KindParameters_Mod
!@sum  This module provides convenient names for standard KIND parameters
!@+    used for declaring integers and floating point numbers.
!@+    This provides greater configurability than the common "real*8" 
!@+    approach, and is supported by the Fortran standard.
!@+    Note that integers should generally be declared _without_
!@+    a kind except where the distinction is important.  This 
!@+    reduces complications when using external libraries.
!@auth T. Clune

  implicit none
  private

  integer, parameter, public :: DP = kind(1.d0)
  integer, parameter, public :: SP = kind(1.0)

  integer, parameter, public :: DEFAULT_INT = kind(1)
  integer, parameter, public :: LONG_INT = selected_int_kind(10)

end module KindParameters_Mod
