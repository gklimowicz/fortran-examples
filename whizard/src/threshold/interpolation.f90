! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'

module interpolation
  use kinds, only: default

  implicit none
  private
  save


  public :: interpolate_linear, strictly_monotonous

  interface interpolate_linear
    module  procedure interpolate_linear_1D_complex_array, &
      interpolate_linear_1D_complex_scalar, &
      interpolate_linear_1D_real_array, &
      interpolate_linear_1D_real_scalar, &
      interpolate_linear_2D_complex_array, &
      interpolate_linear_2D_complex_scalar, &
      interpolate_linear_2D_real_array, &
      interpolate_linear_2D_real_scalar, &
      interpolate_linear_3D_complex_array, &
      interpolate_linear_3D_complex_scalar, &
      interpolate_linear_3D_real_array, &
      interpolate_linear_3D_real_scalar
  end interface

  interface strictly_monotonous
    module procedure monotonous
  end interface strictly_monotonous

  interface find_nearest_left
    !!! recursive bisection is slower
    module procedure find_nearest_left_loop
  end interface find_nearest_left


  interface
    pure module subroutine interpolate_linear_1D_complex_scalar (xa, ya, x, y)
      real(default), dimension(:), intent(in) :: xa
      complex(default), dimension(:), intent(in) :: ya
      real(default), intent(in) :: x
      complex(default), intent(out) :: y
    end subroutine interpolate_linear_1D_complex_scalar
    pure module subroutine interpolate_linear_2D_complex_scalar (x1a, x2a, ya, x1, x2, y)
      real(default), dimension(:), intent(in) :: x1a
      real(default), dimension(:), intent(in) :: x2a
      complex(default), dimension(:,:), intent(in) :: ya
      real(default), intent(in) :: x1
      real(default), intent(in) :: x2
      complex(default), intent(out) :: y
    end subroutine interpolate_linear_2D_complex_scalar
    pure module subroutine interpolate_linear_3D_complex_scalar &
         (x1a, x2a, x3a, ya, x1, x2, x3, y)
      real(default), dimension(:), intent(in) :: x1a
      real(default), dimension(:), intent(in) :: x2a
      real(default), dimension(:), intent(in) :: x3a
      complex(default), dimension(:,:,:), intent(in) :: ya
      real(default), intent(in) :: x1
      real(default), intent(in) :: x2
      real(default), intent(in) :: x3
      complex(default), intent(out) :: y
    end subroutine interpolate_linear_3D_complex_scalar
    pure module subroutine find_nearest_left_loop (xa, x, ixl)
      real(default), dimension(:), intent(in) :: xa
      real(default), intent(in) :: x
      integer, intent(out) :: ixl
    end subroutine find_nearest_left_loop
    pure module function monotonous (xa) result (flag)
      real(default), dimension(:), intent(in) :: xa
      logical :: flag
    end function monotonous
    pure module subroutine interpolate_linear_1D_complex_array (xa, ya, x, y)
      real(default), dimension(:), intent(in) :: xa
      complex(default), dimension(:,:), intent(in) :: ya
      real(default), intent(in) :: x
      complex(default), dimension(size(ya(1,:))), intent(out) :: y
    end subroutine interpolate_linear_1D_complex_array
    pure module subroutine interpolate_linear_1D_real_array (xa, ya, x, y)
      real(default), dimension(:), intent(in) :: xa
      real(default), dimension(:,:), intent(in) :: ya
      real(default), intent(in) :: x
      real(default), dimension(:), intent(out) :: y
    end subroutine interpolate_linear_1D_real_array
    pure module subroutine interpolate_linear_1D_real_scalar (xa, ya, x, y)
      real(default), dimension(:), intent(in) :: xa
      real(default), dimension(:), intent(in) :: ya
      real(default), intent(in) :: x
      real(default), intent(out) :: y
      complex(default), dimension(size(ya)) :: ya_c
    end subroutine interpolate_linear_1D_real_scalar
    pure module subroutine interpolate_linear_2D_complex_array (x1a, x2a, ya, x1, x2, y)
      real(default), dimension(:), intent(in) :: x1a
      real(default), dimension(:), intent(in) :: x2a
      complex(default), dimension(:,:,:), intent(in) :: ya
      real(default), intent(in) :: x1
      real(default), intent(in) :: x2
      complex(default), dimension(size(ya(1,1,:))), intent(out) :: y
    end subroutine interpolate_linear_2D_complex_array
    pure module subroutine interpolate_linear_2D_real_array (x1a, x2a, ya, x1, x2, y)
      real(default), dimension(:), intent(in) :: x1a
      real(default), dimension(:), intent(in) :: x2a
      real(default), dimension(:,:,:), intent(in) :: ya
      real(default), intent(in) :: x1
      real(default), intent(in) :: x2
      real(default), dimension(:), intent(out) :: y
    end subroutine interpolate_linear_2D_real_array
    pure module subroutine interpolate_linear_2D_real_scalar (x1a, x2a, ya, x1, x2, y)
      real(default), dimension(:), intent(in) :: x1a
      real(default), dimension(:), intent(in) :: x2a
      real(default), dimension(:,:), intent(in) :: ya
      real(default), intent(in) :: x1
      real(default), intent(in) :: x2
      real(default), intent(out) :: y
    end subroutine interpolate_linear_2D_real_scalar
    pure module subroutine interpolate_linear_3D_complex_array &
         (x1a, x2a, x3a, ya, x1, x2, x3, y)
      real(default), dimension(:), intent(in) :: x1a
      real(default), dimension(:), intent(in) :: x2a
      real(default), dimension(:), intent(in) :: x3a
      complex(default), dimension(:,:,:,:), intent(in) :: ya
      real(default), intent(in) :: x1
      real(default), intent(in) :: x2
      real(default), intent(in) :: x3
    complex(default), dimension(size(ya(1,1,1,:))), intent(out) :: y
    end subroutine interpolate_linear_3D_complex_array
    pure module subroutine interpolate_linear_3D_real_array (x1a, x2a, x3a, ya, x1, x2, x3, y)
      real(default), dimension(:), intent(in) :: x1a
      real(default), dimension(:), intent(in) :: x2a
      real(default), dimension(:), intent(in) :: x3a
      real(default), dimension(:,:,:,:), intent(in) :: ya
      real(default), intent(in) :: x1
      real(default), intent(in) :: x2
      real(default), intent(in) :: x3
      real(default), dimension(:), intent(out) :: y
    end subroutine interpolate_linear_3D_real_array
    pure module subroutine interpolate_linear_3D_real_scalar (x1a, x2a, x3a, ya, x1, x2, x3, y)
      real(default), dimension(:), intent(in) :: x1a
      real(default), dimension(:), intent(in) :: x2a
      real(default), dimension(:), intent(in) :: x3a
      real(default), dimension(:,:,:), intent(in) :: ya
      real(default), intent(in) :: x1
      real(default), intent(in) :: x2
      real(default), intent(in) :: x3
      real(default), intent(out) :: y
    end subroutine interpolate_linear_3D_real_scalar
  end interface

end module interpolation
