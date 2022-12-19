! coordinates.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module coordinates
  use kinds
  use constants, only: PI
  use specfun, only: gamma
  implicit none
  private
  public :: spherical_to_cartesian_2, &
       spherical_to_cartesian, spherical_to_cartesian_j
  public :: cartesian_to_spherical_2, &
       cartesian_to_spherical, cartesian_to_spherical_j
  public :: spherical_cos_to_cartesian_2, &
       spherical_cos_to_cartesian, spherical_cos_to_cartesian_j
  public :: cartesian_to_spherical_cos_2, &
       cartesian_to_spherical_cos, cartesian_to_spherical_cos_j
  public :: surface
contains
  pure subroutine spherical_to_cartesian_2 (r, phi, theta, x, jacobian)
    real(kind=default), intent(in) :: r, phi
    real(kind=default), dimension(:), intent(in) :: theta
    real(kind=default), dimension(:), intent(out), optional :: x
    real(kind=default), intent(out), optional :: jacobian
    real(kind=default), dimension(size(theta)) :: cos_theta
    real(kind=default), dimension(size(theta)+1) :: product_sin_theta
    integer :: n, i
    n = size (theta) + 2
    cos_theta = cos (theta)
    product_sin_theta(n-1) = 1.0_default
    do i = n - 2, 1, -1
       product_sin_theta(i) = &
            product_sin_theta(i+1) * sqrt (1 - cos_theta(i)**2)
    end do
    if (present (x)) then
       x(1) = r * product_sin_theta(1) * sin (phi)
       x(2) = r * product_sin_theta(1) * cos (phi)
       x(3:) = r * product_sin_theta(2:n-1) * cos_theta
    end if
    if (present (jacobian)) then
       jacobian = r**(n-1) * product (product_sin_theta)
    end if
  end subroutine spherical_to_cartesian_2
  pure function spherical_to_cartesian (r, phi, theta) result (x)
    real(kind=default), intent(in) :: r, phi
    real(kind=default), dimension(:), intent(in) :: theta
    real(kind=default), dimension(size(theta)+2) :: x
    call spherical_to_cartesian_2 (r, phi, theta, x = x)
  end function spherical_to_cartesian
  pure function spherical_to_cartesian_j (r, phi, theta) &
       result (jacobian)
    real(kind=default), intent(in) :: r, phi
    real(kind=default), dimension(:), intent(in) :: theta
    real(kind=default) :: jacobian
    call spherical_to_cartesian_2 (r, phi, theta, jacobian = jacobian)
  end function spherical_to_cartesian_j
  pure subroutine cartesian_to_spherical_2 (x, r, phi, theta, jacobian)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default), intent(out), optional :: r, phi
    real(kind=default), dimension(:), intent(out), optional :: theta
    real(kind=default), intent(out), optional :: jacobian
    real(kind=default) :: local_r
    real(kind=default), dimension(size(x)-2) :: cos_theta
    real(kind=default), dimension(size(x)-1) :: product_sin_theta
    integer :: n, i
    n = size (x)
    local_r = sqrt (dot_product (x, x))
    if (local_r == 0) then
      if (present (r)) then
        r = 0
      end if 
      if (present (phi)) then
        phi = 0
      end if
      if (present (theta)) then
        theta = 0
      end if
      if (present (jacobian)) then
        jacobian = 1
      end if
    else          
      product_sin_theta(n-1) = 1
      do i = n, 3, -1
         if (product_sin_theta(i-1) == 0) then
           cos_theta(i-2) = 0
         else     
           cos_theta(i-2) = x(i) / product_sin_theta(i-1) / local_r
         end if
         product_sin_theta(i-2) = &
              product_sin_theta(i-1) * sqrt (1 - cos_theta(i-2)**2)
      end do
      if (present (r)) then
         r = local_r
      end if
      if (present (phi)) then
         !  Set phi = 0 for vanishing vector
         if (x(1) == 0 .and. x(2)==0) then
          phi = 0
         else     
            phi = atan2 (x(1), x(2))
         end if 
      end if
      if (present (theta)) then
         theta = acos (cos_theta)
      end if
      if (present (jacobian)) then
         jacobian = local_r**(1-n) / product (product_sin_theta)
      end if
    end if
  end subroutine cartesian_to_spherical_2
  pure subroutine cartesian_to_spherical (x, r, phi, theta)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default), intent(out) :: r, phi
    real(kind=default), dimension(:), intent(out) :: theta
    call cartesian_to_spherical_2 (x, r, phi, theta)
  end subroutine cartesian_to_spherical
  pure function cartesian_to_spherical_j (x) result (jacobian)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default) :: jacobian
    call cartesian_to_spherical_2 (x, jacobian = jacobian)
  end function cartesian_to_spherical_j
  pure subroutine spherical_cos_to_cartesian_2 (r, phi, cos_theta, x, jacobian)
    real(kind=default), intent(in) :: r, phi
    real(kind=default), dimension(:), intent(in) :: cos_theta
    real(kind=default), dimension(:), intent(out), optional :: x
    real(kind=default), intent(out), optional :: jacobian
    real(kind=default), dimension(size(cos_theta)+1) :: product_sin_theta
    integer :: n, i
    n = size (cos_theta) + 2
    product_sin_theta(n-1) = 1.0_default
    do i = n - 2, 1, -1
       product_sin_theta(i) = &
            product_sin_theta(i+1) * sqrt (1 - cos_theta(i)**2)
    end do
    if (present (x)) then
       x(1) = r * product_sin_theta(1) * sin (phi)
       x(2) = r * product_sin_theta(1) * cos (phi)
       x(3:) = r * product_sin_theta(2:n-1) * cos_theta
    end if
    if (present (jacobian)) then
       jacobian = r**(n-1) * product (product_sin_theta(2:))
    end if
  end subroutine spherical_cos_to_cartesian_2
  pure function spherical_cos_to_cartesian (r, phi, theta) result (x)
    real(kind=default), intent(in) :: r, phi
    real(kind=default), dimension(:), intent(in) :: theta
    real(kind=default), dimension(size(theta)+2) :: x
    call spherical_cos_to_cartesian_2 (r, phi, theta, x = x)
  end function spherical_cos_to_cartesian
  pure function spherical_cos_to_cartesian_j (r, phi, theta) &
       result (jacobian)
    real(kind=default), intent(in) :: r, phi
    real(kind=default), dimension(:), intent(in) :: theta
    real(kind=default) :: jacobian
    call spherical_cos_to_cartesian_2 (r, phi, theta, jacobian = jacobian)
  end function spherical_cos_to_cartesian_j
  pure subroutine cartesian_to_spherical_cos_2 (x, r, phi, cos_theta, jacobian)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default), intent(out), optional :: r, phi
    real(kind=default), dimension(:), intent(out), optional :: cos_theta
    real(kind=default), intent(out), optional :: jacobian
    real(kind=default) :: local_r
    real(kind=default), dimension(size(x)-2) :: local_cos_theta
    real(kind=default), dimension(size(x)-1) :: product_sin_theta
    integer :: n, i
    n = size (x)
    local_r = sqrt (dot_product (x, x))
    if (local_r == 0) then
      if (present (r)) then
        r = 0
      end if 
      if (present (phi)) then
        phi = 0
      end if
      if (present (cos_theta)) then
        cos_theta = 0
      end if
      if (present (jacobian)) then
        jacobian = 1
      end if
    else          
      product_sin_theta(n-1) = 1
      do i = n, 3, -1
         if (product_sin_theta(i-1) == 0) then
           local_cos_theta(i-2) = 0
         else     
           local_cos_theta(i-2) = x(i) / product_sin_theta(i-1) / local_r
         end if
         product_sin_theta(i-2) = &
              product_sin_theta(i-1) * sqrt (1 - local_cos_theta(i-2)**2)
      end do
      if (present (r)) then
         r = local_r
      end if
      if (present (phi)) then
         !  Set phi = 0 for vanishing vector
         if (x(1) == 0 .and. x(2)==0) then
          phi = 0
         else     
            phi = atan2 (x(1), x(2))
         end if 
      end if
      if (present (cos_theta)) then
         cos_theta = local_cos_theta
      end if
      if (present (jacobian)) then
         jacobian = local_r**(1-n) / product (product_sin_theta(2:))
      end if
    end if
  end subroutine cartesian_to_spherical_cos_2
  pure subroutine cartesian_to_spherical_cos (x, r, phi, cos_theta)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default), intent(out) :: r, phi
    real(kind=default), dimension(:), intent(out), optional :: cos_theta
    call cartesian_to_spherical_cos_2 (x, r, phi, cos_theta)
  end subroutine cartesian_to_spherical_cos
  pure function cartesian_to_spherical_cos_j (x) result (jacobian)
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default) :: jacobian
    call cartesian_to_spherical_cos_2 (x, jacobian = jacobian)
  end function cartesian_to_spherical_cos_j
  pure function surface (n) result (vol)
    integer, intent(in) :: n
    real(kind=default) :: vol
    real(kind=default) :: n_by_2
    n_by_2 = 0.5_default * n 
    vol = 2 * PI**n_by_2 / gamma (n_by_2)
  end function surface
end module coordinates
