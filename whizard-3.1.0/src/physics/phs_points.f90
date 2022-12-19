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

module phs_points

  use kinds, only: default
  use lorentz, only: vector4_t
  use lorentz, only: lorentz_transformation_t
  use lorentz, only: sum

  implicit none
  private

  public :: phs_point_t
  public :: assignment(=)
  public :: size
  public :: operator(==)
  public :: operator(*)
  public :: sum

  type :: phs_point_t
     private
     type(vector4_t), dimension(:), allocatable :: p
  contains
    procedure :: write => phs_point_write
    procedure :: get => phs_point_get
    procedure :: select => phs_point_select
    procedure :: get_msq => phs_point_get_msq
    procedure :: get_x => phs_point_get_x
  end type phs_point_t


  interface assignment(=)
     module procedure phs_point_from_n
     module procedure phs_point_from_vector4
     module procedure vector4_from_phs_point
  end interface
  interface size
     module procedure phs_point_size
  end interface size
  interface operator(==)
     module procedure phs_point_eq
  end interface operator(==)
  interface operator(*)
     module procedure prod_LT_phs_point
  end interface operator(*)
  interface sum
     module procedure phs_point_sum
     module procedure phs_point_sum_iarray
  end interface sum

  interface
    module subroutine phs_point_write (phs_point, unit, show_mass, testflag, &
        check_conservation, ultra, n_in)
      class(phs_point_t), intent(in) :: phs_point
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_mass
      logical, intent(in), optional :: testflag, ultra
      logical, intent(in), optional :: check_conservation
      integer, intent(in), optional :: n_in
    end subroutine phs_point_write
    pure module subroutine phs_point_from_n (phs_point, n_particles)
      type(phs_point_t), intent(out) :: phs_point
      integer, intent(in) :: n_particles
    end subroutine phs_point_from_n
    pure module subroutine phs_point_from_vector4 (phs_point, p)
      type(phs_point_t), intent(out) :: phs_point
      type(vector4_t), dimension(:), intent(in) :: p
    end subroutine phs_point_from_vector4
    pure module subroutine vector4_from_phs_point (p, phs_point)
      class(phs_point_t), intent(in) :: phs_point
      type(vector4_t), dimension(:), allocatable, intent(out) :: p
    end subroutine vector4_from_phs_point
    pure module function phs_point_size (phs_point) result (s)
      class(phs_point_t), intent(in) :: phs_point
      integer :: s
    end function phs_point_size
    elemental module function phs_point_eq &
         (phs_point_1, phs_point_2) result (flag)
      class(phs_point_t), intent(in) :: phs_point_1, phs_point_2
      logical :: flag
    end function phs_point_eq
    pure module function phs_point_get (phs_point) result (p)
      class(phs_point_t), intent(in) :: phs_point
      type(vector4_t), dimension(:), allocatable :: p
    end function phs_point_get
    elemental module function phs_point_select (phs_point, i) result (p)
      class(phs_point_t), intent(in) :: phs_point
      integer, intent(in) :: i
      type(vector4_t) :: p
    end function phs_point_select
    pure module function phs_point_get_msq (phs_point, iarray) result (msq)
      class(phs_point_t), intent(in) :: phs_point
      integer, dimension(:), intent(in) :: iarray
      real(default) :: msq
    end function phs_point_get_msq
    elemental module function prod_LT_phs_point (L, phs_point) result (phs_point_LT)
      type(lorentz_transformation_t), intent(in) :: L
      type(phs_point_t), intent(in) :: phs_point
      type(phs_point_t) :: phs_point_LT
    end function prod_LT_phs_point
    pure module function phs_point_sum (phs_point, mask) result (p)
      class(phs_point_t), intent(in) :: phs_point
      logical, dimension(:), intent(in), optional :: mask
      type(vector4_t) :: p
    end function phs_point_sum
    pure module function phs_point_sum_iarray (phs_point, iarray) result (p)
      class(phs_point_t), intent(in) :: phs_point
      integer, dimension(:), intent(in) :: iarray
      type(vector4_t) :: p
    end function phs_point_sum_iarray
    pure module function phs_point_get_x (phs_point, E_beam) result (x)
      class(phs_point_t), intent(in) :: phs_point
      real(default), dimension(2) :: x
      real(default), intent(in) :: E_beam
    end function phs_point_get_x
  end interface

end module phs_points
