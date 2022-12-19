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

module helicities

  implicit none
  private

  public :: helicity_t
  public :: helicity

  type :: helicity_t
     private
     logical :: defined = .false.
     integer :: h1, h2
   contains
     generic :: init => helicity_init_empty, helicity_init_same, helicity_init_different
     procedure, private :: helicity_init_empty
     procedure, private :: helicity_init_same
     procedure, private :: helicity_init_different
     procedure :: undefine => helicity_undefine
     procedure :: diagonalize => helicity_diagonalize
     procedure :: flip => helicity_flip
     procedure :: get_indices => helicity_get_indices
     procedure :: write => helicity_write
     procedure :: write_raw => helicity_write_raw
     procedure :: read_raw => helicity_read_raw
     procedure :: is_defined => helicity_is_defined
     procedure :: is_diagonal => helicity_is_diagonal
     procedure :: to_pair => helicity_to_pair
     generic :: operator(.match.) => helicity_match
     generic :: operator(.dmatch.) => helicity_match_diagonal
     generic :: operator(==) => helicity_eq
     generic :: operator(/=) => helicity_neq
     procedure, private ::  helicity_match
     procedure, private ::  helicity_match_diagonal
     procedure, private ::  helicity_eq
     procedure, private ::  helicity_neq
     generic :: operator(.merge.) => merge_helicities
     procedure, private ::  merge_helicities
  end type helicity_t


  interface helicity
     module procedure helicity0, helicity1, helicity2
  end interface helicity


  interface
    pure module function helicity0 () result (hel)
      type(helicity_t) :: hel
    end function helicity0
    elemental module function helicity1 (h) result (hel)
      type(helicity_t) :: hel
      integer, intent(in) :: h
    end function helicity1
    elemental module function helicity2 (h2, h1) result (hel)
      type(helicity_t) :: hel
      integer, intent(in) :: h1, h2
    end function helicity2
    elemental module subroutine helicity_init_empty (hel)
      class(helicity_t), intent(inout) :: hel
    end subroutine helicity_init_empty
    elemental module subroutine helicity_init_same (hel, h)
      class(helicity_t), intent(inout) :: hel
      integer, intent(in) :: h
    end subroutine helicity_init_same
    elemental module subroutine helicity_init_different (hel, h2, h1)
      class(helicity_t), intent(inout) :: hel
      integer, intent(in) :: h1, h2
    end subroutine helicity_init_different
    elemental module subroutine helicity_undefine (hel)
      class(helicity_t), intent(inout) :: hel
    end subroutine helicity_undefine
    elemental module subroutine helicity_diagonalize (hel)
      class(helicity_t), intent(inout) :: hel
    end subroutine helicity_diagonalize
    elemental module subroutine helicity_flip (hel)
      class(helicity_t), intent(inout) :: hel
    end subroutine helicity_flip
    module subroutine helicity_get_indices (hel, h1, h2)
      class(helicity_t), intent(in) :: hel
      integer, intent(out) :: h1, h2
    end subroutine helicity_get_indices
    module subroutine helicity_write (hel, unit)
      class(helicity_t), intent(in) :: hel
      integer, intent(in), optional :: unit
    end subroutine helicity_write
    module subroutine helicity_write_raw (hel, u)
      class(helicity_t), intent(in) :: hel
      integer, intent(in) :: u
    end subroutine helicity_write_raw
    module subroutine helicity_read_raw (hel, u, iostat)
      class(helicity_t), intent(out) :: hel
      integer, intent(in) :: u
      integer, intent(out), optional :: iostat
    end subroutine helicity_read_raw
    elemental module function helicity_is_defined (hel) result (defined)
      logical :: defined
      class(helicity_t), intent(in) :: hel
    end function helicity_is_defined
    elemental module function helicity_is_diagonal (hel) result (diagonal)
      logical :: diagonal
      class(helicity_t), intent(in) :: hel
    end function helicity_is_diagonal
    pure module function helicity_to_pair (hel) result (h)
      integer, dimension(2) :: h
      class(helicity_t), intent(in) :: hel
    end function helicity_to_pair
    elemental module function helicity_match (hel1, hel2) result (eq)
      logical :: eq
      class(helicity_t), intent(in) :: hel1, hel2
    end function helicity_match
    elemental module function helicity_match_diagonal (hel1, hel2) result (eq)
      logical :: eq
      class(helicity_t), intent(in) :: hel1, hel2
    end function helicity_match_diagonal
    elemental module function helicity_eq (hel1, hel2) result (eq)
      logical :: eq
      class(helicity_t), intent(in) :: hel1, hel2
    end function helicity_eq
    elemental module function helicity_neq (hel1, hel2) result (neq)
      logical :: neq
      class(helicity_t), intent(in) :: hel1, hel2
    end function helicity_neq
    elemental module function merge_helicities (hel1, hel2) result (hel)
      type(helicity_t) :: hel
      class(helicity_t), intent(in) :: hel1, hel2
    end function merge_helicities
  end interface

end module helicities
