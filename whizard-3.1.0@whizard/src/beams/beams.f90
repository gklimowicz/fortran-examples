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

module beams

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use model_data
  use flavors
  use quantum_numbers
  use state_matrices
  use interactions
  use polarizations
  use beam_structures

  implicit none
  private

  public :: beam_data_t
  public :: beam_t
  public :: beam_init
  public :: beam_final
  public :: beam_write
  public :: assignment(=)
  public :: interaction_set_source_link_beam
  public :: beam_get_int_ptr
  public :: beam_set_momenta

  type :: beam_data_t
     logical :: initialized = .false.
     integer :: n = 0
     type(flavor_t), dimension(:), allocatable :: flv
     real(default), dimension(:), allocatable :: mass
     type(pmatrix_t), dimension(:), allocatable :: pmatrix
     logical :: lab_is_cm = .true.
     type(vector4_t), dimension(:), allocatable :: p_cm
     type(vector4_t), dimension(:), allocatable :: p
     type(lorentz_transformation_t), allocatable  :: L_cm_to_lab
     real(default) :: sqrts = 0
     character(32) :: md5sum = ""
   contains
     procedure :: final => beam_data_final
     procedure :: write => beam_data_write
     procedure :: are_valid => beam_data_are_valid
     procedure :: check_scattering => beam_data_check_scattering
     procedure :: get_n_in => beam_data_get_n_in
     procedure :: get_flavor => beam_data_get_flavor
     procedure :: get_energy => beam_data_get_energy
     procedure :: get_sqrts => beam_data_get_sqrts
     procedure :: get_polarization => beam_data_get_polarization
     procedure :: get_helicity_state_matrix => beam_data_get_helicity_state_matrix
     procedure :: is_initialized => beam_data_is_initialized
     procedure :: get_md5sum => beam_data_get_md5sum
     procedure :: init_structure => beam_data_init_structure
     procedure :: init_sqrts => beam_data_init_sqrts
     procedure :: init_momenta => beam_data_init_momenta
     procedure :: compute_md5sum => beam_data_compute_md5sum
     procedure :: init_decay => beam_data_init_decay
  end type beam_data_t

  type :: beam_t
     private
     type(interaction_t) :: int
  end type beam_t


  interface assignment(=)
     module procedure beam_assign
  end interface


  interface
    module subroutine beam_data_final (beam_data)
      class(beam_data_t), intent(inout) :: beam_data
    end subroutine beam_data_final
    module subroutine beam_data_write (beam_data, unit, verbose, write_md5sum)
      class(beam_data_t), intent(in) :: beam_data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, write_md5sum
    end subroutine beam_data_write
    module function beam_data_are_valid (beam_data) result (flag)
      class(beam_data_t), intent(in) :: beam_data
      logical :: flag
    end function beam_data_are_valid
    module subroutine beam_data_check_scattering (beam_data, sqrts)
      class(beam_data_t), intent(in) :: beam_data
      real(default), intent(in), optional :: sqrts
    end subroutine beam_data_check_scattering
    module function beam_data_get_n_in (beam_data) result (n_in)
      class(beam_data_t), intent(in) :: beam_data
      integer :: n_in
    end function beam_data_get_n_in
    module function beam_data_get_flavor (beam_data) result (flv)
      class(beam_data_t), intent(in) :: beam_data
      type(flavor_t), dimension(:), allocatable :: flv
    end function beam_data_get_flavor
    module function beam_data_get_energy (beam_data) result (e)
      class(beam_data_t), intent(in) :: beam_data
      real(default), dimension(:), allocatable :: e
    end function beam_data_get_energy
    module function beam_data_get_sqrts (beam_data) result (sqrts)
      class(beam_data_t), intent(in) :: beam_data
      real(default) :: sqrts
    end function beam_data_get_sqrts
    module function beam_data_get_polarization (beam_data) result (pol)
      class(beam_data_t), intent(in) :: beam_data
      real(default), dimension(beam_data%n) :: pol
    end function beam_data_get_polarization
    module function beam_data_get_helicity_state_matrix &
         (beam_data) result (state_hel)
      type(state_matrix_t) :: state_hel
      class(beam_data_t), intent(in) :: beam_data
    end function beam_data_get_helicity_state_matrix
    module function beam_data_is_initialized (beam_data) result (initialized)
      logical :: initialized
      class(beam_data_t), intent(in) :: beam_data
    end function beam_data_is_initialized
    module function beam_data_get_md5sum &
         (beam_data, sqrts) result (md5sum_beams)
      class(beam_data_t), intent(in) :: beam_data
      real(default), intent(in) :: sqrts
      character(32) :: md5sum_beams
    end function beam_data_get_md5sum
    module subroutine beam_data_init_structure &
         (beam_data, structure, sqrts, model, decay_rest_frame)
      class(beam_data_t), intent(out) :: beam_data
      type(beam_structure_t), intent(in) :: structure
      real(default), intent(in) :: sqrts
      class(model_data_t), intent(in), target :: model
      logical, intent(in), optional :: decay_rest_frame
    end subroutine beam_data_init_structure
    module subroutine beam_data_init_sqrts &
         (beam_data, sqrts, flv, smatrix, pol_f)
      class(beam_data_t), intent(out) :: beam_data
      real(default), intent(in) :: sqrts
      type(flavor_t), dimension(:), intent(in) :: flv
      type(smatrix_t), dimension(:), intent(in), optional :: smatrix
      real(default), dimension(:), intent(in), optional :: pol_f
    end subroutine beam_data_init_sqrts
    module subroutine beam_data_init_momenta &
         (beam_data, p3, flv, smatrix, pol_f)
      class(beam_data_t), intent(out) :: beam_data
      type(vector3_t), dimension(:), intent(in) :: p3
      type(flavor_t), dimension(:), intent(in) :: flv
      type(smatrix_t), dimension(:), intent(in), optional :: smatrix
      real(default), dimension(:), intent(in), optional :: pol_f
    end subroutine beam_data_init_momenta
    module subroutine beam_data_compute_md5sum (beam_data)
      class(beam_data_t), intent(inout) :: beam_data
      integer :: unit
    end subroutine beam_data_compute_md5sum
    module subroutine beam_data_init_decay &
         (beam_data, flv, smatrix, pol_f, rest_frame)
      class(beam_data_t), intent(out) :: beam_data
      type(flavor_t), dimension(1), intent(in) :: flv
      type(smatrix_t), dimension(1), intent(in), optional :: smatrix
      real(default), dimension(:), intent(in), optional :: pol_f
      logical, intent(in), optional :: rest_frame
    end subroutine beam_data_init_decay
    module subroutine beam_init (beam, beam_data)
      type(beam_t), intent(out) :: beam
      type(beam_data_t), intent(in), target :: beam_data
    end subroutine beam_init
    module subroutine beam_final (beam)
      type(beam_t), intent(inout) :: beam
    end subroutine beam_final
    module subroutine beam_write &
         (beam, unit, verbose, show_momentum_sum, show_mass, col_verbose)
      type(beam_t), intent(in) :: beam
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, show_momentum_sum, show_mass
      logical, intent(in), optional :: col_verbose
    end subroutine beam_write
    module subroutine beam_assign (beam_out, beam_in)
      type(beam_t), intent(out) :: beam_out
      type(beam_t), intent(in) :: beam_in
    end subroutine beam_assign
    module subroutine interaction_set_source_link_beam (int, i, beam1, i1)
      type(interaction_t), intent(inout) :: int
      type(beam_t), intent(in), target :: beam1
      integer, intent(in) :: i, i1
    end subroutine interaction_set_source_link_beam
    module function beam_get_int_ptr (beam) result (int)
      type(interaction_t), pointer :: int
      type(beam_t), intent(in), target :: beam
    end function beam_get_int_ptr
    module subroutine beam_set_momenta (beam, p)
      type(beam_t), intent(inout) :: beam
      type(vector4_t), dimension(:), intent(in) :: p
    end subroutine beam_set_momenta
  end interface

end module beams
