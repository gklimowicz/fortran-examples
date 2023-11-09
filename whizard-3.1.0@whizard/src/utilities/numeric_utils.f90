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

module numeric_utils

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: assert
  public:: assert_equal
  interface assert_equal
     module procedure assert_equal_integer, assert_equal_integers, &
            assert_equal_real, assert_equal_reals, &
            assert_equal_complex, assert_equal_complexs
  end interface

  public :: nearly_equal
  public:: vanishes
  interface vanishes
     module procedure vanishes_real, vanishes_complex
  end interface
  public :: expanded_amp2
  public :: abs2
  public:: remove_array_element
  interface remove_array_element
     module procedure remove_array_element_logical
  end interface
  public :: remove_duplicates_from_int_array
  public :: extend_integer_array
  public :: crop_integer_array
  public :: log_prec
  public :: split_array
  public :: d1mach
  public :: dqk41
  public :: dqk61
  public :: gauss_kronrod
  public :: pacify

  integer, parameter, public :: GAUSS_KRONROD_41 = 41, GAUSS_KRONROD_61 = 61

  type :: int_workspace_t
     private
     integer :: limit
     integer :: size = 0
     integer :: nrmax = 1
     integer :: i = 1
     integer :: maximum_level = 0
     real(default), dimension(:), allocatable :: alist
     real(default), dimension(:), allocatable :: blist
     real(default), dimension(:), allocatable :: rlist
     real(default), dimension(:), allocatable :: elist
     integer, dimension(:), allocatable :: order
     integer, dimension(:), allocatable :: level
   contains
     procedure :: init => int_workspace_init
     procedure :: set_initial => int_workspace_set_initial
     procedure :: update => int_workspace_update
     procedure :: sort => int_workspace_sort
  end type int_workspace_t


  interface nearly_equal
     module procedure nearly_equal_real
     module procedure nearly_equal_complex
  end interface nearly_equal

  interface split_array
     module procedure split_integer_array
     module procedure split_real_array
  end interface
  abstract interface
     function g_func (x) result (f)
       import default
       real(default), intent(in) :: x
       real(default) :: f
     end function g_func
  end interface
  interface pacify
     module procedure pacify_real_default
     module procedure pacify_complex_default
  end interface pacify


  interface
    module subroutine int_workspace_init (work, a, b, limit)
      class(int_workspace_t), intent(out) :: work
      real(default), intent(in) :: a, b
      integer, intent(in) :: limit
    end subroutine int_workspace_init
    module subroutine int_workspace_set_initial (work, res, err)
      class(int_workspace_t), intent(inout) :: work
      real(default), intent(in) :: res, err
    end subroutine int_workspace_set_initial
    module subroutine int_workspace_update (work, a1, b1, area1, error1, &
         a2, b2, area2, error2)
      class(int_workspace_t), intent(inout) :: work
      real(default), intent(in) :: a1, b1, area1, error1, &
           a2, b2, area2, error2
    end subroutine int_workspace_update
    module subroutine int_workspace_sort (work)
      class(int_workspace_t), intent(inout) :: work
    end subroutine int_workspace_sort
    module subroutine assert (unit, ok, description, exit_on_fail)
      integer, intent(in) :: unit
      logical, intent(in) :: ok
      character(*), intent(in), optional :: description
      logical, intent(in), optional :: exit_on_fail
    end subroutine assert
    module subroutine assert_equal_integer (unit, lhs, rhs, description, exit_on_fail)
      integer, intent(in) :: unit
      integer, intent(in) :: lhs, rhs
      character(*), intent(in), optional :: description
      logical, intent(in), optional :: exit_on_fail
    end subroutine assert_equal_integer
    module subroutine assert_equal_integers (unit, lhs, rhs, description, exit_on_fail)
      integer, intent(in) :: unit
      integer, dimension(:), intent(in) :: lhs, rhs
      character(*), intent(in), optional :: description
      logical, intent(in), optional :: exit_on_fail
    end subroutine assert_equal_integers
    module subroutine assert_equal_real (unit, lhs, rhs, description, &
                                  abs_smallness, rel_smallness, exit_on_fail)
      integer, intent(in) :: unit
      real(default), intent(in) :: lhs, rhs
      character(*), intent(in), optional :: description
      real(default), intent(in), optional :: abs_smallness, rel_smallness
      logical, intent(in), optional :: exit_on_fail
    end subroutine assert_equal_real
    module subroutine assert_equal_reals (unit, lhs, rhs, description, &
                                  abs_smallness, rel_smallness, exit_on_fail)
      integer, intent(in) :: unit
      real(default), dimension(:), intent(in) :: lhs, rhs
      character(*), intent(in), optional :: description
      real(default), intent(in), optional :: abs_smallness, rel_smallness
      logical, intent(in), optional :: exit_on_fail
    end subroutine assert_equal_reals
    module subroutine assert_equal_complex (unit, lhs, rhs, description, &
                                  abs_smallness, rel_smallness, exit_on_fail)
      integer, intent(in) :: unit
      complex(default), intent(in) :: lhs, rhs
      character(*), intent(in), optional :: description
      real(default), intent(in), optional :: abs_smallness, rel_smallness
      logical, intent(in), optional :: exit_on_fail
    end subroutine assert_equal_complex
    module subroutine assert_equal_complexs (unit, lhs, rhs, description, &
                                  abs_smallness, rel_smallness, exit_on_fail)
      integer, intent(in) :: unit
      complex(default), dimension(:), intent(in) :: lhs, rhs
      character(*), intent(in), optional :: description
      real(default), intent(in), optional :: abs_smallness, rel_smallness
      logical, intent(in), optional :: exit_on_fail
    end subroutine assert_equal_complexs
    elemental module function nearly_equal_real &
         (a, b, abs_smallness, rel_smallness) result (r)
      logical :: r
      real(default), intent(in) :: a, b
      real(default), intent(in), optional :: abs_smallness, rel_smallness
    end function nearly_equal_real
    elemental module function nearly_equal_complex &
         (a, b, abs_smallness, rel_smallness) result (r)
      logical :: r
      complex(default), intent(in) :: a, b
      real(default), intent(in), optional :: abs_smallness, rel_smallness
    end function nearly_equal_complex
    elemental module function vanishes_real &
         (x, abs_smallness, rel_smallness) result (r)
      logical :: r
      real(default), intent(in) :: x
      real(default), intent(in), optional :: abs_smallness, rel_smallness
    end function vanishes_real
    elemental module function vanishes_complex &
         (x, abs_smallness, rel_smallness) result (r)
      logical :: r
      complex(default), intent(in) :: x
      real(default), intent(in), optional :: abs_smallness, rel_smallness
    end function vanishes_complex
    pure module function expanded_amp2 (amp_tree, amp_blob) result (amp2)
      real(default) :: amp2
      complex(default), dimension(:), intent(in) :: amp_tree, amp_blob
    end function expanded_amp2
    elemental module function abs2 (c) result (c2)
      real(default) :: c2
      complex(default), intent(in) :: c
    end function abs2
    module function remove_array_element_logical &
         (array, index) result (array_reduced)
      logical, intent(in), dimension(:) :: array
      integer, intent(in) :: index
      logical, dimension(:), allocatable :: array_reduced
    end function remove_array_element_logical
    module function remove_duplicates_from_int_array &
         (array) result (array_unique)
      integer, intent(in), dimension(:) :: array
      integer, dimension(:), allocatable :: array_unique
    end function remove_duplicates_from_int_array
    module subroutine extend_integer_array (list, incr, initial_value)
      integer, intent(inout), dimension(:), allocatable :: list
      integer, intent(in) :: incr
      integer, intent(in), optional :: initial_value
    end subroutine extend_integer_array
    module subroutine crop_integer_array (list, i_crop)
      integer, intent(inout), dimension(:), allocatable :: list
      integer, intent(in) :: i_crop
    end subroutine crop_integer_array
    module function log_prec (x, xb) result (lx)
      real(default), intent(in) :: x, xb
      real(default) :: lx
    end function log_prec
    module subroutine split_integer_array (list1, list2)
      integer, intent(inout), dimension(:), allocatable :: list1, list2
      integer, dimension(:), allocatable :: list_store
    end subroutine split_integer_array
    module subroutine split_real_array (list1, list2)
      real(default), intent(inout), dimension(:), allocatable :: list1, list2
      real(default), dimension(:), allocatable :: list_store
    end subroutine split_real_array
    module function d1mach (i) result (d1)
      integer, intent(in) :: i
      real(default) :: d1
    end function d1mach
    module subroutine dqk41 (f, a, b, result, abserr, resabs, resasc)
      procedure(g_func) :: f
      real(default), intent(in) :: a, b
      real(default), intent(out) :: result, abserr, resabs, resasc
    end subroutine dqk41
    module subroutine dqk61 (f, a, b, result, abserr, resabs, resasc)
      procedure(g_func) :: f
      real(default), intent(in) :: a, b
      real(default), intent(out) :: result, abserr, resabs, resasc
    end subroutine dqk61
   module subroutine gauss_kronrod &
        (type, f, a, b, limit, result, abserr, epsabs, epsrel)
     integer, intent(in) :: type
     procedure(g_func) :: f
     real(default), intent(in) :: a, b
     real(default), intent(in) :: epsabs, epsrel
     real(default), intent(out) :: result, abserr
     integer, intent(in) :: limit
   end subroutine gauss_kronrod
    elemental module subroutine pacify_real_default (x, tolerance)
      real(default), intent(inout) :: x
      real(default), intent(in) :: tolerance
    end subroutine pacify_real_default

    elemental module subroutine pacify_complex_default (x, tolerance)
      complex(default), intent(inout) :: x
      real(default), intent(in) :: tolerance
    end subroutine pacify_complex_default
  end interface

end module numeric_utils
