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

submodule (interpolation) interpolation_s

  implicit none



contains

  pure module subroutine interpolate_linear_1D_complex_scalar (xa, ya, x, y)
    real(default), dimension(:), intent(in) :: xa
    complex(default), dimension(:), intent(in) :: ya
    real(default), intent(in) :: x
    complex(default), intent(out) :: y
    integer :: ixl
    real(default) :: t
    y = 0.0_default
    !!! don't check this at runtime:
    ! if ( .not.monotonous(xa) ) return
    if ( out_of_range(xa, x) ) return
    ixl = 0
    call find_nearest_left (xa, x, ixl)
    t = ( x - xa(ixl) ) / ( xa(ixl+1) - xa(ixl) )
    y = (1.-t)*ya(ixl) + t*ya(ixl+1)
  end subroutine interpolate_linear_1D_complex_scalar

  pure module subroutine interpolate_linear_2D_complex_scalar (x1a, x2a, ya, x1, x2, y)
    real(default), dimension(:), intent(in) :: x1a
    real(default), dimension(:), intent(in) :: x2a
    complex(default), dimension(:,:), intent(in) :: ya
    real(default), intent(in) :: x1
    real(default), intent(in) :: x2
    complex(default), intent(out) :: y
    integer :: ix1l, ix2l
    real(default) :: t, u
    y = 0.0_default
    !!! don't check this at runtime:
    ! if ( (.not.monotonous(x1a)) .or. (.not.monotonous(x2a)) ) return
    if ( out_of_range(x1a, x1) .or. out_of_range(x2a, x2) ) return
    ix1l = 0
    call find_nearest_left (x1a, x1, ix1l)
    ix2l = 0
    call find_nearest_left (x2a, x2, ix2l)
    t = ( x1 - x1a(ix1l) ) / ( x1a(ix1l+1) - x1a(ix1l) )
    u = ( x2 - x2a(ix2l) ) / ( x2a(ix2l+1) - x2a(ix2l) )
    y =  (1.-t)*(1.-u)*ya(ix1l  ,ix2l  ) &
        +    t *(1.-u)*ya(ix1l+1,ix2l  ) &
        +    t *    u *ya(ix1l+1,ix2l+1) &
        +(1.-t)*    u *ya(ix1l  ,ix2l+1)
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
    integer :: ix1l, ix2l, ix3l
    real(default) :: t, u, v
    y = 0.0_default
    !!! don't check this at runtime:
    ! if ( (.not.monotonous(x1a)) .or. (.not.monotonous(x2a)) ) return
    if ( out_of_range(x1a, x1) .or. out_of_range(x2a, x2) .or. out_of_range(x3a, x3) ) return
    ix1l = 0
    call find_nearest_left (x1a, x1, ix1l)
    ix2l = 0
    call find_nearest_left (x2a, x2, ix2l)
    ix3l = 0
    call find_nearest_left (x3a, x3, ix3l)
    t = ( x1 - x1a(ix1l) ) / ( x1a(ix1l+1) - x1a(ix1l) )
    u = ( x2 - x2a(ix2l) ) / ( x2a(ix2l+1) - x2a(ix2l) )
    v = ( x3 - x3a(ix3l) ) / ( x3a(ix3l+1) - x3a(ix3l) )
    y =  (1.-t)*(1.-u)*(1.-v)*ya(ix1l  ,ix2l  ,ix3l  ) &
        +(1.-t)*(1.-u)*    v *ya(ix1l  ,ix2l  ,ix3l+1) &
        +(1.-t)*    u *(1.-v)*ya(ix1l  ,ix2l+1,ix3l  ) &
        +(1.-t)*    u *    v *ya(ix1l  ,ix2l+1,ix3l+1) &
        +    t *(1.-u)*(1.-v)*ya(ix1l+1,ix2l  ,ix3l  ) &
        +    t *(1.-u)*    v *ya(ix1l+1,ix2l  ,ix3l+1) &
        +    t *    u *(1.-v)*ya(ix1l+1,ix2l+1,ix3l  ) &
        +    t *    u *    v *ya(ix1l+1,ix2l+1,ix3l+1)
  end subroutine interpolate_linear_3D_complex_scalar

  pure module subroutine find_nearest_left_loop (xa, x, ixl)
    real(default), dimension(:), intent(in) :: xa
    real(default), intent(in) :: x
    integer, intent(out) :: ixl
    integer :: ixm, ixr
    ixl = 1
    ixr = size(xa)
    do
      if ( ixr-ixl <= 1 ) return
      ixm = (ixr+ixl) / 2
      if ( x < xa(ixm) ) then
        ixr = ixm
      else
        ixl = ixm
      end if
    end do
  end subroutine find_nearest_left_loop

  pure recursive subroutine find_nearest_left_rec (xa, x, ixl)
    real(default), dimension(:), intent(in) :: xa
    real(default), intent(in) :: x
    integer, intent(inout) :: ixl
    integer :: nx, bs
    real(default), dimension(:), allocatable :: xa_new
    nx = size(xa)
    bs = nx/2 + 1
    if ( nx < 3 ) then
      ixl = ixl + bs - 1
      return
    else
      if ( x < xa(bs) ) then
        allocate( xa_new(1:bs) )
        xa_new = xa(1:bs)
      else
        ixl = ixl + bs - 1
        allocate( xa_new(bs:nx) )
        xa_new = xa(bs:nx)
      end if
      call find_nearest_left_rec (xa_new, x, ixl)
      deallocate( xa_new )
    end if
  end subroutine find_nearest_left_rec

  pure module function monotonous (xa) result (flag)
    real(default), dimension(:), intent(in) :: xa
    integer :: ix
    logical :: flag
    flag = .false.
    do ix = 1, size(xa)-1
      flag = ( xa(ix) < xa(ix+1) )
      if ( .not. flag ) return
    end do
  end function monotonous

  pure function out_of_range (xa, x) result (flag)
    real(default), dimension(:), intent(in) :: xa
    real(default), intent(in) :: x
    logical :: flag
    flag = ( x < xa(1) .or. x > xa(size(xa)) )
  end function out_of_range

  pure module subroutine interpolate_linear_1D_complex_array (xa, ya, x, y)
    real(default), dimension(:), intent(in) :: xa
    complex(default), dimension(:,:), intent(in) :: ya
    real(default), intent(in) :: x
    complex(default), dimension(size(ya(1,:))), intent(out) :: y
    integer :: iy
    do iy=1, size(y)
      call interpolate_linear_1D_complex_scalar (xa, ya(:,iy), x, y(iy))
    end do
  end subroutine interpolate_linear_1D_complex_array

  pure module subroutine interpolate_linear_1D_real_array (xa, ya, x, y)
    real(default), dimension(:), intent(in) :: xa
    real(default), dimension(:,:), intent(in) :: ya
    real(default), intent(in) :: x
    real(default), dimension(:), intent(out) :: y
    complex(default), dimension(size(ya(1,:))) :: y_c
    call interpolate_linear_1D_complex_array (xa, cmplx(ya,kind=default), x, y_c)
    y = real(y_c,kind=default)
  end subroutine interpolate_linear_1D_real_array

  pure module subroutine interpolate_linear_1D_real_scalar (xa, ya, x, y)
    real(default), dimension(:), intent(in) :: xa
    real(default), dimension(:), intent(in) :: ya
    real(default), intent(in) :: x
    real(default), intent(out) :: y
    complex(default), dimension(size(ya)) :: ya_c
    complex(default) :: y_c
    ya_c = cmplx(ya,kind=default)
    call interpolate_linear_1D_complex_scalar (xa, ya_c, x, y_c)
    y = real(y_c,kind=default)
  end subroutine interpolate_linear_1D_real_scalar

  pure module subroutine interpolate_linear_2D_complex_array (x1a, x2a, ya, x1, x2, y)
    real(default), dimension(:), intent(in) :: x1a
    real(default), dimension(:), intent(in) :: x2a
    complex(default), dimension(:,:,:), intent(in) :: ya
    real(default), intent(in) :: x1
    real(default), intent(in) :: x2
    complex(default), dimension(size(ya(1,1,:))), intent(out) :: y
    integer :: iy
    do iy=1, size(y)
      call interpolate_linear_2D_complex_scalar (x1a, x2a, ya(:,:,iy), x1, x2, y(iy))
    end do
  end subroutine interpolate_linear_2D_complex_array

  pure module subroutine interpolate_linear_2D_real_array (x1a, x2a, ya, x1, x2, y)
    real(default), dimension(:), intent(in) :: x1a
    real(default), dimension(:), intent(in) :: x2a
    real(default), dimension(:,:,:), intent(in) :: ya
    real(default), intent(in) :: x1
    real(default), intent(in) :: x2
    real(default), dimension(:), intent(out) :: y
    complex(default), dimension(size(ya(1,1,:))) :: y_c
    call interpolate_linear_2D_complex_array (x1a, x2a, cmplx(ya,kind=default), x1, x2, y_c)
    y = real(y_c,kind=default)
  end subroutine interpolate_linear_2D_real_array

  pure module subroutine interpolate_linear_2D_real_scalar (x1a, x2a, ya, x1, x2, y)
    real(default), dimension(:), intent(in) :: x1a
    real(default), dimension(:), intent(in) :: x2a
    real(default), dimension(:,:), intent(in) :: ya
    real(default), intent(in) :: x1
    real(default), intent(in) :: x2
    real(default), intent(out) :: y
    complex(default), dimension(size(ya(:,1)),size(ya(1,:))) :: ya_c
    complex(default) :: y_c
    ya_c = reshape (ya_c, shape(ya))
    ya_c = cmplx(ya,kind=default)
    call interpolate_linear_2D_complex_scalar (x1a, x2a, ya_c, x1, x2, y_c)
    y = real(y_c,kind=default)
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
    integer :: iy
    do iy=1, size(y)
      call interpolate_linear_3D_complex_scalar &
                            (x1a, x2a, x3a, ya(:,:,:,iy), x1, x2, x3, y(iy))
    end do
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
    complex(default), dimension(size(ya(1,1,1,:))) :: y_c
    call interpolate_linear_3D_complex_array &
                       (x1a, x2a, x3a, cmplx(ya,kind=default), x1, x2, x3, y_c)
    y = real(y_c,kind=default)
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
    complex(default), dimension(size(ya(:,1,1)),size(ya(1,:,1)),size(ya(1,1,:))) :: ya_c
    complex(default) :: y_c
    ya_c = cmplx(ya,kind=default)
    call interpolate_linear_3D_complex_scalar (x1a, x2a, x3a, ya_c, x1, x2, x3, y_c)
    y = real(y_c,kind=default)
  end subroutine interpolate_linear_3D_real_scalar


end submodule interpolation_s

