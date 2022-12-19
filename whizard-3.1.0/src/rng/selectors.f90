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

module selectors

  use kinds, only: default
  use rng_base

  implicit none
  private

  public :: selector_t

  type :: selector_t
     integer :: offset = 0
     integer, dimension(:), allocatable :: map
     real(default), dimension(:), allocatable :: weight
     real(default), dimension(:), allocatable :: acc
   contains
     procedure :: write => selector_write
     procedure :: init => selector_init
     procedure :: select => selector_select
     procedure :: generate => selector_generate
     procedure :: get_weight => selector_get_weight
  end type selector_t


  interface
    module subroutine selector_write (object, unit, testflag)
      class(selector_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine selector_write
    module subroutine selector_init (selector, weight, negative_weights, offset)
      class(selector_t), intent(out) :: selector
      real(default), dimension(:), intent(in) :: weight
      logical, intent(in), optional :: negative_weights
      integer, intent(in), optional :: offset
    end subroutine selector_init
    module function selector_select (selector, x) result (n)
      class(selector_t), intent(in) :: selector
      real(default), intent(in) :: x
      integer :: n
    end function selector_select
    module subroutine selector_generate (selector, rng, n)
      class(selector_t), intent(in) :: selector
      class(rng_t), intent(inout) :: rng
      integer, intent(out) :: n
    end subroutine selector_generate
    module function selector_get_weight (selector, n) result (weight)
      class(selector_t), intent(in) :: selector
      integer, intent(in) :: n
      real(default) :: weight
    end function selector_get_weight
  end interface

end module selectors
