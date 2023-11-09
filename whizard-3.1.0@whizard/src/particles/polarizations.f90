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

module polarizations

  use kinds, only: default
  use physics_defs, only: SCALAR, SPINOR, VECTOR, VECTORSPINOR, TENSOR
  use flavors
  use quantum_numbers
  use state_matrices
  use bloch_vectors

  implicit none
  private

  public :: polarization_t
  public :: combine_polarization_states
  public :: polarization_iterator_t
  public :: smatrix_t
  public :: pmatrix_t

  type :: polarization_t
     private
     integer :: spin_type = SCALAR
     integer :: multiplicity = 1
     integer :: chirality = 0
     logical :: anti = .false.
     type(bloch_vector_t) :: bv
   contains
     generic, private :: init => polarization_init, polarization_init_flv
     procedure, private :: polarization_init
     procedure, private :: polarization_init_flv
     generic :: init_generic => &
          polarization_init_generic, &
          polarization_init_generic_flv
     procedure, private :: polarization_init_generic
     procedure, private :: polarization_init_generic_flv
     procedure :: write => polarization_write
     procedure :: write_raw => polarization_write_raw
     procedure :: read_raw => polarization_read_raw
     procedure :: is_polarized => polarization_is_polarized
     procedure :: is_diagonal => polarization_is_diagonal
     procedure :: init_state_matrix => polarization_init_state_matrix
     procedure :: to_state => polarization_to_state_matrix
     procedure :: init_unpolarized => polarization_init_unpolarized
     procedure :: init_circular => polarization_init_circular
     procedure :: init_transversal => polarization_init_transversal
     procedure :: init_axis => polarization_init_axis
     procedure :: init_angles => polarization_init_angles
     procedure :: init_longitudinal => polarization_init_longitudinal
     procedure :: init_diagonal => polarization_init_diagonal
     procedure :: get_axis => polarization_get_axis
     procedure :: to_angles => polarization_to_angles
     procedure :: init_pmatrix => polarization_init_pmatrix
  end type polarization_t

  type :: polarization_iterator_t
     private
     type(polarization_t), pointer :: pol => null ()
     logical :: polarized = .false.
     integer :: h1 = 0
     integer :: h2 = 0
     integer :: i = 0
     integer :: j = 0
     complex(default), dimension(:,:), allocatable :: r
     complex(default) :: value = 1._default
     real(default) :: tolerance = -1._default
     logical :: valid = .false.
   contains
     procedure :: write => polarization_iterator_write
     procedure :: init => polarization_iterator_init
     procedure :: advance => polarization_iterator_advance
     procedure :: is_valid => polarization_iterator_is_valid
     procedure :: get_value => polarization_iterator_get_value
     procedure :: get_quantum_numbers => polarization_iterator_get_quantum_numbers
  end type polarization_iterator_t

  type :: smatrix_t
     private
     integer :: dim = 0
     integer :: n_entry = 0
     integer, dimension(:,:), allocatable :: index
     complex(default), dimension(:), allocatable :: value
   contains
     procedure :: write => smatrix_write
     procedure :: init => smatrix_init
     procedure :: set_entry => smatrix_set_entry
     procedure :: exists => smatrix_exists
  end type smatrix_t

  type, extends (smatrix_t) :: pmatrix_t
     private
     integer :: spin_type = 0
     integer :: multiplicity = 0
     logical :: massive = .true.
     integer :: chirality = 0
     real(default) :: degree = 1
     logical :: pure = .false.
   contains
     procedure :: write => pmatrix_write
     generic :: assignment(=) => pmatrix_assign_from_smatrix
     procedure, private :: pmatrix_assign_from_smatrix
     procedure :: normalize => pmatrix_normalize
     procedure :: is_polarized => pmatrix_is_polarized
     procedure :: is_diagonal => pmatrix_is_diagonal
     procedure :: get_simple_pol => pmatrix_get_simple_pol
  end type pmatrix_t




  interface
    module subroutine polarization_init (pol, spin_type, multiplicity, &
         anti, left_handed, right_handed)
      class(polarization_t), intent(out) :: pol
      integer, intent(in) :: spin_type
      integer, intent(in) :: multiplicity
      logical, intent(in) :: anti
      logical, intent(in) :: left_handed
      logical, intent(in) :: right_handed
    end subroutine polarization_init
    module subroutine polarization_init_flv (pol, flv)
      class(polarization_t), intent(out) :: pol
      type(flavor_t), intent(in) :: flv
    end subroutine polarization_init_flv
    module subroutine polarization_init_generic (pol, spin_type, multiplicity, &
         anti, left_handed, right_handed)
      class(polarization_t), intent(out) :: pol
      integer, intent(in) :: spin_type
      integer, intent(in) :: multiplicity
      logical, intent(in) :: anti
      logical, intent(in) :: left_handed
      logical, intent(in) :: right_handed
    end subroutine polarization_init_generic
    module subroutine polarization_init_generic_flv (pol, flv)
      class(polarization_t), intent(out) :: pol
      type(flavor_t), intent(in) :: flv
    end subroutine polarization_init_generic_flv
    module subroutine polarization_write (pol, unit, state_matrix, all_states, tolerance)
      class(polarization_t), intent(in) :: pol
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: state_matrix, all_states
      real(default), intent(in), optional :: tolerance
    end subroutine polarization_write
    module subroutine polarization_write_raw (pol, u)
      class(polarization_t), intent(in) :: pol
      integer, intent(in) :: u
    end subroutine polarization_write_raw
    module subroutine polarization_read_raw (pol, u, iostat)
      class(polarization_t), intent(out) :: pol
      integer, intent(in) :: u
      integer, intent(out), optional :: iostat
    end subroutine polarization_read_raw
    module function polarization_is_polarized (pol) result (polarized)
      class(polarization_t), intent(in) :: pol
      logical :: polarized
    end function polarization_is_polarized
    module function polarization_is_diagonal (pol) result (diagonal)
      class(polarization_t), intent(in) :: pol
      logical :: diagonal
    end function polarization_is_diagonal
    module subroutine polarization_init_state_matrix (pol, state)
      class(polarization_t), intent(out) :: pol
      type(state_matrix_t), intent(in), target :: state
    end subroutine polarization_init_state_matrix
    module subroutine polarization_to_state_matrix (pol, state, all_states, tolerance)
      class(polarization_t), intent(in), target :: pol
      type(state_matrix_t), intent(out) :: state
      logical, intent(in), optional :: all_states
      real(default), intent(in), optional :: tolerance
    end subroutine polarization_to_state_matrix
    module subroutine polarization_init_unpolarized (pol, flv)
      class(polarization_t), intent(out) :: pol
      type(flavor_t), intent(in) :: flv
    end subroutine polarization_init_unpolarized
    module subroutine polarization_init_circular (pol, flv, f)
      class(polarization_t), intent(out) :: pol
      type(flavor_t), intent(in) :: flv
      real(default), intent(in) :: f
    end subroutine polarization_init_circular
    module subroutine polarization_init_transversal (pol, flv, phi, f)
      class(polarization_t), intent(out) :: pol
      type(flavor_t), intent(in) :: flv
      real(default), intent(in) :: phi, f
    end subroutine polarization_init_transversal
    module subroutine polarization_init_axis (pol, flv, alpha)
      class(polarization_t), intent(out) :: pol
      type(flavor_t), intent(in) :: flv
      real(default), dimension(3), intent(in) :: alpha
    end subroutine polarization_init_axis
    module subroutine polarization_init_angles (pol, flv, r, theta, phi)
      class(polarization_t), intent(out) :: pol
      type(flavor_t), intent(in) :: flv
      real(default), intent(in) :: r, theta, phi
    end subroutine polarization_init_angles
    module subroutine polarization_init_longitudinal (pol, flv, f)
      class(polarization_t), intent(out) :: pol
      type(flavor_t), intent(in) :: flv
      real(default), intent(in) :: f
    end subroutine polarization_init_longitudinal
    module subroutine polarization_init_diagonal (pol, flv, rd)
      class(polarization_t), intent(out) :: pol
      type(flavor_t), intent(in) :: flv
      real(default), dimension(:), intent(in) :: rd
    end subroutine polarization_init_diagonal
    module subroutine combine_polarization_states (pol, state)
      type(polarization_t), dimension(:), intent(in), target :: pol
      type(state_matrix_t), intent(out) :: state
    end subroutine combine_polarization_states
    module function polarization_get_axis (pol) result (alpha)
      class(polarization_t), intent(in), target :: pol
      real(default), dimension(3) :: alpha
    end function polarization_get_axis
    module subroutine polarization_to_angles (pol, r, theta, phi)
      class(polarization_t), intent(in) :: pol
      real(default), intent(out) :: r, theta, phi
    end subroutine polarization_to_angles
    module subroutine polarization_iterator_write (it, unit)
      class(polarization_iterator_t), intent(in) :: it
      integer, intent(in), optional :: unit
    end subroutine polarization_iterator_write
    module subroutine polarization_iterator_init (it, pol, all_states, tolerance)
      class(polarization_iterator_t), intent(out) :: it
      type(polarization_t), intent(in), target :: pol
      logical, intent(in), optional :: all_states
      real(default), intent(in), optional :: tolerance
    end subroutine polarization_iterator_init
    recursive module subroutine polarization_iterator_advance (it)
      class(polarization_iterator_t), intent(inout) :: it
    end subroutine polarization_iterator_advance
    module function polarization_iterator_is_valid (it) result (is_valid)
      logical :: is_valid
      class(polarization_iterator_t), intent(in) :: it
    end function polarization_iterator_is_valid
    module function polarization_iterator_get_value (it) result (value)
      complex(default) :: value
      class(polarization_iterator_t), intent(in) :: it
    end function polarization_iterator_get_value
    module function polarization_iterator_get_quantum_numbers (it) result (qn)
      class(polarization_iterator_t), intent(in) :: it
      type(quantum_numbers_t) :: qn
    end function polarization_iterator_get_quantum_numbers
    module subroutine smatrix_write (object, unit, indent)
      class(smatrix_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
    end subroutine smatrix_write
    module subroutine smatrix_init (smatrix, dim, n_entry)
      class(smatrix_t), intent(out) :: smatrix
      integer, intent(in) :: dim
      integer, intent(in) :: n_entry
    end subroutine smatrix_init
    module subroutine smatrix_set_entry (smatrix, i, index, value)
      class(smatrix_t), intent(inout) :: smatrix
      integer, intent(in) :: i
      integer, dimension(:), intent(in) :: index
      complex(default), intent(in) :: value
    end subroutine smatrix_set_entry
    elemental module function smatrix_exists (smatrix) result (exist)
      logical :: exist
      class(smatrix_t), intent(in) :: smatrix
    end function smatrix_exists
    module subroutine pmatrix_write (object, unit, indent)
      class(pmatrix_t), intent(in) :: object
      integer, intent(in), optional :: unit, indent
    end subroutine pmatrix_write
    module subroutine pmatrix_assign_from_smatrix (pmatrix, smatrix)
      class(pmatrix_t), intent(out) :: pmatrix
      type(smatrix_t), intent(in) :: smatrix
    end subroutine pmatrix_assign_from_smatrix
    module subroutine pmatrix_normalize (pmatrix, flv, degree, tolerance)
      class(pmatrix_t), intent(inout) :: pmatrix
      type(flavor_t), intent(in) :: flv
      real(default), intent(in), optional :: degree
      real(default), intent(in), optional :: tolerance
    end subroutine pmatrix_normalize
    elemental module function pmatrix_is_polarized (pmatrix) result (flag)
      class(pmatrix_t), intent(in) :: pmatrix
      logical :: flag
    end function pmatrix_is_polarized
    elemental module function pmatrix_is_diagonal (pmatrix) result (flag)
      class(pmatrix_t), intent(in) :: pmatrix
      logical :: flag
    end function pmatrix_is_diagonal
    elemental module function pmatrix_get_simple_pol (pmatrix) result (pol)
      class(pmatrix_t), intent(in) :: pmatrix
      real(default) :: pol
    end function pmatrix_get_simple_pol
    module subroutine polarization_init_pmatrix (pol, pmatrix)
      class(polarization_t), intent(out) :: pol
      type(pmatrix_t), intent(in) :: pmatrix
    end subroutine polarization_init_pmatrix
  end interface

end module polarizations
