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

module sf_beam_events

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use file_registries
  use pdg_arrays
  use model_data
  use flavors
  use quantum_numbers
  use state_matrices
  use polarizations
  use sf_base

  implicit none
  private

  public :: beam_file_registry
  public :: beam_events_data_t
  public :: beam_events_t

  type, extends(sf_data_t) :: beam_events_data_t
     private
     type(flavor_t), dimension(2) :: flv_in
     type(string_t) :: dir
     type(string_t) :: file
     type(string_t) :: fqn
     integer :: unit = 0
     logical :: warn_eof = .true.
   contains
     procedure :: init => beam_events_data_init
     procedure :: is_generator => beam_events_data_is_generator
     procedure :: get_n_par => beam_events_data_get_n_par
     procedure :: get_pdg_out => beam_events_data_get_pdg_out
     procedure :: allocate_sf_int => beam_events_data_allocate_sf_int
     procedure :: write => beam_events_data_write
     procedure :: open => beam_events_data_open
     procedure :: close => beam_events_data_close
     procedure :: get_beam_file => beam_events_data_get_beam_file
  end type beam_events_data_t

  type, extends (sf_int_t) :: beam_events_t
     type(beam_events_data_t), pointer :: data => null ()
     integer :: count = 0
   contains
     procedure :: type_string => beam_events_type_string
     procedure :: write => beam_events_write
     procedure :: init => beam_events_init
     procedure :: final => sf_beam_events_final
     procedure :: is_generator => beam_events_is_generator
     procedure :: generate_free => beam_events_generate_free
     procedure :: complete_kinematics => beam_events_complete_kinematics
     procedure :: inverse_kinematics => beam_events_inverse_kinematics
     procedure :: apply => beam_events_apply
  end type beam_events_t


  type(file_registry_t), save :: beam_file_registry


  interface
    module subroutine beam_events_data_init &
         (data, model, pdg_in, dir, file, warn_eof)
      class(beam_events_data_t), intent(out) :: data
      class(model_data_t), intent(in), target :: model
      type(pdg_array_t), dimension(2), intent(in) :: pdg_in
      type(string_t), intent(in) :: dir
      type(string_t), intent(in) :: file
      logical, intent(in), optional :: warn_eof
    end subroutine beam_events_data_init
    module function beam_events_data_is_generator (data) result (flag)
      class(beam_events_data_t), intent(in) :: data
      logical :: flag
    end function beam_events_data_is_generator
    module function beam_events_data_get_n_par (data) result (n)
      class(beam_events_data_t), intent(in) :: data
      integer :: n
    end function beam_events_data_get_n_par
    module subroutine beam_events_data_get_pdg_out (data, pdg_out)
      class(beam_events_data_t), intent(in) :: data
      type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    end subroutine beam_events_data_get_pdg_out
    module subroutine beam_events_data_write (data, unit, verbose)
      class(beam_events_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine beam_events_data_write
    module subroutine beam_events_data_open (data)
      class(beam_events_data_t), intent(inout) :: data
    end subroutine beam_events_data_open
    module subroutine beam_events_data_close (data)
      class(beam_events_data_t), intent(inout) :: data
    end subroutine beam_events_data_close
    module function beam_events_data_get_beam_file (data) result (file)
      class(beam_events_data_t), intent(in) :: data
      type(string_t) :: file
    end function beam_events_data_get_beam_file
    module function beam_events_type_string (object) result (string)
      class(beam_events_t), intent(in) :: object
      type(string_t) :: string
    end function beam_events_type_string
    module subroutine beam_events_write (object, unit, testflag)
      class(beam_events_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine beam_events_write
    module subroutine beam_events_init (sf_int, data)
      class(beam_events_t), intent(out) :: sf_int
      class(sf_data_t), intent(in), target :: data
    end subroutine beam_events_init
    module subroutine sf_beam_events_final (object)
      class(beam_events_t), intent(inout) :: object
    end subroutine sf_beam_events_final
    module function beam_events_is_generator (sf_int) result (flag)
      class(beam_events_t), intent(in) :: sf_int
      logical :: flag
    end function beam_events_is_generator
    recursive module subroutine beam_events_generate_free &
         (sf_int, r, rb,  x_free)
      class(beam_events_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(inout) :: x_free
    end subroutine beam_events_generate_free
    module subroutine beam_events_complete_kinematics &
         (sf_int, x, xb, f, r, rb, map)
      class(beam_events_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), intent(in) :: rb
      logical, intent(in) :: map
    end subroutine beam_events_complete_kinematics
    module subroutine beam_events_inverse_kinematics &
         (sf_int, x, xb, f, r, rb, map, set_momenta)
      class(beam_events_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: r
      real(default), dimension(:), intent(out) :: rb
      logical, intent(in) :: map
      logical, intent(in), optional :: set_momenta
    end subroutine beam_events_inverse_kinematics
    module subroutine beam_events_apply &
         (sf_int, scale, negative_sf, rescale, i_sub)
      class(beam_events_t), intent(inout) :: sf_int
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(in), optional :: rescale
      integer, intent(in), optional :: i_sub
    end subroutine beam_events_apply
  end interface

contains

  subroutine beam_events_data_allocate_sf_int (data, sf_int)
    class(beam_events_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (beam_events_t :: sf_int)
  end subroutine beam_events_data_allocate_sf_int


end module sf_beam_events
