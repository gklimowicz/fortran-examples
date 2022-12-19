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

module iterations

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: iterations_list_t
  public :: iteration_multipliers_t

  type :: iterations_spec_t
     private
     integer :: n_it = 0
     integer :: n_calls = 0
     logical :: custom_adaptation = .false.
     logical :: adapt_grids = .false.
     logical :: adapt_weights = .false.
  end type iterations_spec_t

  type :: iterations_list_t
     private
     integer :: n_pass = 0
     type(iterations_spec_t), dimension(:), allocatable :: pass
   contains
     procedure :: init => iterations_list_init
     procedure :: clear => iterations_list_clear
     procedure :: write => iterations_list_write
     procedure :: to_string => iterations_list_to_string
     procedure :: get_n_pass => iterations_list_get_n_pass
     procedure :: get_n_calls => iterations_list_get_n_calls
     procedure :: set_n_calls => iterations_list_set_n_calls
     procedure :: adapt_grids => iterations_list_adapt_grids
     procedure :: adapt_weights => iterations_list_adapt_weights
     procedure :: get_n_it => iterations_list_get_n_it
  end type iterations_list_t

  type :: iteration_multipliers_t
    real(default) :: mult_real = 1._default
    real(default) :: mult_virt = 1._default
    real(default) :: mult_dglap = 1._default
    real(default) :: mult_threshold = 1._default
    integer, dimension(:), allocatable :: n_calls0
  end type iteration_multipliers_t


  interface
    module subroutine iterations_list_init &
         (it_list, n_it, n_calls, adapt, adapt_code, adapt_grids, adapt_weights)
      class(iterations_list_t), intent(inout) :: it_list
      integer, dimension(:), intent(in) :: n_it, n_calls
      logical, dimension(:), intent(in), optional :: adapt
      type(string_t), dimension(:), intent(in), optional :: adapt_code
      logical, dimension(:), intent(in), optional :: adapt_grids, adapt_weights
    end subroutine iterations_list_init
    module subroutine iterations_list_clear (it_list)
      class(iterations_list_t), intent(inout) :: it_list
    end subroutine iterations_list_clear
    module subroutine iterations_list_write (it_list, unit)
      class(iterations_list_t), intent(in) :: it_list
      integer, intent(in), optional :: unit
    end subroutine iterations_list_write
    module function iterations_list_to_string (it_list) result (buffer)
      class(iterations_list_t), intent(in) :: it_list
      type(string_t) :: buffer
    end function iterations_list_to_string
    module function iterations_list_get_n_pass (it_list) result (n_pass)
      class(iterations_list_t), intent(in) :: it_list
      integer :: n_pass
    end function iterations_list_get_n_pass
    module function iterations_list_get_n_calls (it_list, pass) result (n_calls)
      class(iterations_list_t), intent(in) :: it_list
      integer :: n_calls
      integer, intent(in) :: pass
    end function iterations_list_get_n_calls
    module subroutine iterations_list_set_n_calls (it_list, pass, n_calls)
      class(iterations_list_t), intent(inout) :: it_list
      integer, intent(in) :: pass, n_calls
    end subroutine iterations_list_set_n_calls
    module function iterations_list_adapt_grids (it_list, pass) result (flag)
      logical :: flag
      class(iterations_list_t), intent(in) :: it_list
      integer, intent(in) :: pass
    end function iterations_list_adapt_grids
    module function iterations_list_adapt_weights (it_list, pass) result (flag)
      logical :: flag
      class(iterations_list_t), intent(in) :: it_list
      integer, intent(in) :: pass
    end function iterations_list_adapt_weights
    module function iterations_list_get_n_it (it_list, pass) result (n_it)
      class(iterations_list_t), intent(in) :: it_list
      integer :: n_it
      integer, intent(in) :: pass
    end function iterations_list_get_n_it
  end interface

end module iterations
