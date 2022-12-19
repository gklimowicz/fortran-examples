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

submodule (phs_points) phs_points_s

  use lorentz, only: vector4_null
  use lorentz, only: vector4_write_set
  use lorentz, only: operator(==)
  use lorentz, only: operator(*)
  use lorentz, only: operator(**)

  implicit none

contains

  module subroutine phs_point_write (phs_point, unit, show_mass, testflag, &
      check_conservation, ultra, n_in)
    class(phs_point_t), intent(in) :: phs_point
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_mass
    logical, intent(in), optional :: testflag, ultra
    logical, intent(in), optional :: check_conservation
    integer, intent(in), optional :: n_in
    if (allocated (phs_point%p)) then
       call vector4_write_set (phs_point%p, &
            unit = unit, &
            show_mass = show_mass, &
            testflag = testflag, &
            check_conservation = check_conservation, &
            ultra = ultra, &
            n_in = n_in)
    end if
  end subroutine phs_point_write

  pure module subroutine phs_point_from_n (phs_point, n_particles)
    type(phs_point_t), intent(out) :: phs_point
    integer, intent(in) :: n_particles
    allocate (phs_point%p (n_particles), source = vector4_null)
  end subroutine phs_point_from_n

  pure module subroutine phs_point_from_vector4 (phs_point, p)
    type(phs_point_t), intent(out) :: phs_point
    type(vector4_t), dimension(:), intent(in) :: p
    phs_point%p = p
  end subroutine phs_point_from_vector4

  pure module subroutine vector4_from_phs_point (p, phs_point)
    class(phs_point_t), intent(in) :: phs_point
    type(vector4_t), dimension(:), allocatable, intent(out) :: p
    if (allocated (phs_point%p))  p = phs_point%p
  end subroutine vector4_from_phs_point

  pure module function phs_point_size (phs_point) result (s)
    class(phs_point_t), intent(in) :: phs_point
    integer :: s
    if (allocated (phs_point%p)) then
       s = size (phs_point%p)
    else
       s = 0
    end if
  end function phs_point_size

  elemental module function phs_point_eq &
       (phs_point_1, phs_point_2) result (flag)
    class(phs_point_t), intent(in) :: phs_point_1, phs_point_2
    logical :: flag
    if (allocated (phs_point_1%p) .and. (allocated (phs_point_2%p))) then
       flag = all (phs_point_1%p == phs_point_2%p)
    else
       flag = .false.
    end if
  end function phs_point_eq

  pure module function phs_point_get (phs_point) result (p)
    class(phs_point_t), intent(in) :: phs_point
    type(vector4_t), dimension(:), allocatable :: p
    if (allocated (phs_point%p)) then
       p = phs_point%p
    else
       allocate (p (0))
    end if
  end function phs_point_get

  elemental module function phs_point_select (phs_point, i) result (p)
    class(phs_point_t), intent(in) :: phs_point
    integer, intent(in) :: i
    type(vector4_t) :: p
    if (allocated (phs_point%p)) then
       p = phs_point%p(i)
    else
       p = vector4_null
    end if
  end function phs_point_select

  pure module function phs_point_get_msq (phs_point, iarray) result (msq)
    class(phs_point_t), intent(in) :: phs_point
    integer, dimension(:), intent(in) :: iarray
    real(default) :: msq
    if (allocated (phs_point%p)) then
       msq = (sum (phs_point%p(iarray)))**2
    else
       msq = 0
    end if
  end function phs_point_get_msq

  elemental module function prod_LT_phs_point (L, phs_point) result (phs_point_LT)
    type(lorentz_transformation_t), intent(in) :: L
    type(phs_point_t), intent(in) :: phs_point
    type(phs_point_t) :: phs_point_LT
    if (allocated (phs_point%p))  phs_point_LT%p = L * phs_point%p
  end function prod_LT_phs_point

  pure module function phs_point_sum (phs_point, mask) result (p)
    class(phs_point_t), intent(in) :: phs_point
    logical, dimension(:), intent(in), optional :: mask
    type(vector4_t) :: p
    if (allocated (phs_point%p)) then
       if (present (mask)) then
          p = sum (phs_point%p, mask)
       else
          p = sum (phs_point%p)
       end if
    else
       p = vector4_null
    end if
  end function phs_point_sum

  pure module function phs_point_sum_iarray (phs_point, iarray) result (p)
    class(phs_point_t), intent(in) :: phs_point
    integer, dimension(:), intent(in) :: iarray
    type(vector4_t) :: p
    logical, dimension(:), allocatable :: mask
    integer :: i
    allocate (mask (size (phs_point)), source = .false.)
    mask(iarray) = .true.
    p = sum (phs_point, mask)
  end function phs_point_sum_iarray

  pure module function phs_point_get_x (phs_point, E_beam) result (x)
    class(phs_point_t), intent(in) :: phs_point
    real(default), dimension(2) :: x
    real(default), intent(in) :: E_beam
    x = phs_point%p(1:2)%p(0) / E_beam
  end function phs_point_get_x


end submodule phs_points_s

