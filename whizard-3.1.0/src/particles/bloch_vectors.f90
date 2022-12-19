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

module bloch_vectors

  use kinds, only: default
  use physics_defs, only: UNKNOWN

  implicit none
  private

  public :: bloch_vector_t

  type :: bloch_vector_t
     private
     integer :: spin_type = UNKNOWN
     real(default), dimension(:), allocatable :: a
   contains
     procedure :: init_unpolarized => bloch_vector_init_unpolarized
     generic :: init => bloch_vector_init
     procedure, private :: bloch_vector_init
     procedure :: from_array => bloch_vector_from_array
     procedure :: to_array => bloch_vector_to_array
     procedure :: write_raw => bloch_vector_write_raw
     procedure :: read_raw => bloch_vector_read_raw
     procedure :: get_n_states
     procedure :: get_length
     procedure :: hel_index => bv_helicity_index
     procedure :: hel_value => bv_helicity_value
     procedure :: bloch_factor => bv_factor
     procedure :: is_defined => bloch_vector_is_defined
     procedure :: is_polarized => bloch_vector_is_polarized
     procedure :: is_diagonal => bloch_vector_is_diagonal
     procedure :: get_norm => bloch_vector_get_norm
     generic :: init => bloch_vector_init_diagonal
     procedure, private :: bloch_vector_init_diagonal
     generic :: set => bloch_vector_set_diagonal
     procedure, private :: bloch_vector_set_diagonal
     procedure :: init_max_weight => bloch_vector_init_max_weight
     procedure :: init_vector => bloch_vector_init_vector
     procedure :: to_vector => bloch_vector_to_vector
     generic :: init => bloch_vector_init_matrix
     procedure, private :: bloch_vector_init_matrix
     generic :: set => bloch_vector_set_matrix
     procedure, private :: bloch_vector_set_matrix
     procedure :: to_matrix => bloch_vector_to_matrix
end type bloch_vector_t


  interface
    module subroutine bloch_vector_init_unpolarized (pol, spin_type)
      class(bloch_vector_t), intent(out) :: pol
      integer, intent(in) :: spin_type
    end subroutine bloch_vector_init_unpolarized
    module subroutine bloch_vector_init (pol, spin_type)
      class(bloch_vector_t), intent(out) :: pol
      integer, intent(in) :: spin_type
    end subroutine bloch_vector_init
    module subroutine bloch_vector_from_array (pol, a)
      class(bloch_vector_t), intent(inout) :: pol
      real(default), dimension(:), allocatable, intent(in) :: a
    end subroutine bloch_vector_from_array
    module subroutine bloch_vector_to_array (pol, a)
      class(bloch_vector_t), intent(in) :: pol
      real(default), dimension(:), allocatable, intent(out) :: a
    end subroutine bloch_vector_to_array
    module subroutine bloch_vector_write_raw (pol, u)
      class(bloch_vector_t), intent(in) :: pol
      integer, intent(in) :: u
    end subroutine bloch_vector_write_raw
    module subroutine bloch_vector_read_raw (pol, u, iostat)
      class(bloch_vector_t), intent(out) :: pol
      integer, intent(in) :: u
      integer, intent(out) :: iostat
    end subroutine bloch_vector_read_raw
    module function get_n_states (pol) result (n)
      class(bloch_vector_t), intent(in) :: pol
      integer :: n
    end function get_n_states
    module function get_length (pol) result (n)
      class(bloch_vector_t), intent(in) :: pol
      integer :: n
    end function get_length
    module function bv_helicity_index (pol, h) result (i)
      class(bloch_vector_t), intent(in) :: pol
      integer, intent(in) :: h
      integer :: i
    end function bv_helicity_index
    module function bv_helicity_value (pol, i) result (h)
      class(bloch_vector_t), intent(in) :: pol
      integer, intent(in) :: i
      integer :: h
    end function bv_helicity_value
    module function bv_factor (pol) result (f)
      class(bloch_vector_t), intent(in) :: pol
      real(default) :: f
    end function bv_factor
    module function bloch_vector_is_defined (pol) result (flag)
      class(bloch_vector_t), intent(in) :: pol
      logical :: flag
    end function bloch_vector_is_defined
    module function bloch_vector_is_polarized (pol) result (flag)
      class(bloch_vector_t), intent(in) :: pol
      logical :: flag
    end function bloch_vector_is_polarized
    module function bloch_vector_is_diagonal (pol) result (diagonal)
      class(bloch_vector_t), intent(in) :: pol
      logical :: diagonal
    end function bloch_vector_is_diagonal
    module function bloch_vector_get_norm (pol) result (norm)
      class(bloch_vector_t), intent(in) :: pol
      real(default) :: norm
    end function bloch_vector_get_norm
    module subroutine bloch_vector_init_diagonal (pol, spin_type, rd)
      class(bloch_vector_t), intent(out) :: pol
      integer, intent(in) :: spin_type
      real(default), dimension(:), intent(in) :: rd
    end subroutine bloch_vector_init_diagonal
    module subroutine bloch_vector_set_diagonal (pol, rd)
      class(bloch_vector_t), intent(inout) :: pol
      real(default), dimension(:), intent(in) :: rd
    end subroutine bloch_vector_set_diagonal
    module subroutine bloch_vector_init_max_weight (pol, spin_type)
      class(bloch_vector_t), intent(out) :: pol
      integer, intent(in) :: spin_type
    end subroutine bloch_vector_init_max_weight
    module subroutine bloch_vector_init_vector (pol, s, a)
      class(bloch_vector_t), intent(out) :: pol
      integer, intent(in) :: s
      real(default), dimension(3), intent(in) :: a
    end subroutine bloch_vector_init_vector
    module subroutine bloch_vector_to_vector (pol, a)
      class(bloch_vector_t), intent(in) :: pol
      real(default), dimension(3), intent(out) :: a
    end subroutine bloch_vector_to_vector
    module subroutine bloch_vector_init_matrix (pol, spin_type, r)
      class(bloch_vector_t), intent(out) :: pol
      integer, intent(in) :: spin_type
      complex(default), dimension(:,:), intent(in) :: r
    end subroutine bloch_vector_init_matrix
    module subroutine bloch_vector_set_matrix (pol, r)
      class(bloch_vector_t), intent(inout) :: pol
      complex(default), dimension(:,:), intent(in) :: r
    end subroutine bloch_vector_set_matrix
    module subroutine bloch_vector_to_matrix (pol, r, only_max_weight)
      class(bloch_vector_t), intent(in) :: pol
      complex(default), dimension(:,:), intent(out), allocatable :: r
      logical, intent(in), optional :: only_max_weight
    end subroutine bloch_vector_to_matrix
  end interface

end module bloch_vectors
