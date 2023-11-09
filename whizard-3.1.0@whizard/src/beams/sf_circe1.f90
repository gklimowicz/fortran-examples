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

module sf_circe1

  use kinds, only: default
  use kinds, only: double
  use iso_varying_string, string_t => varying_string
  use rng_base
  use pdg_arrays
  use model_data
  use flavors
  use polarizations
  use sf_mappings
  use sf_base
  use circe1, circe1_rng_t => rng_type !NODEP!

  implicit none
  private

  public :: circe1_data_t
  public :: circe1_t

  type, extends (sf_data_t) :: circe1_data_t
     private
     class(model_data_t), pointer :: model => null ()
     type(flavor_t), dimension(2) :: flv_in
     integer, dimension(2) :: pdg_in
     real(default), dimension(2) :: m_in = 0
     logical, dimension(2) :: photon = .false.
     logical :: generate = .false.
     class(rng_factory_t), allocatable :: rng_factory
     real(default) :: sqrts = 0
     real(default) :: eps = 0
     integer :: ver = 0
     integer :: rev = 0
     character(6) :: acc = "?"
     integer :: chat = 0
     logical :: with_radiation = .false.
   contains
       procedure :: init => circe1_data_init
       procedure :: set_generator_mode => circe1_data_set_generator_mode
       procedure :: check => circe1_data_check
       procedure :: write => circe1_data_write
       procedure :: is_generator => circe1_data_is_generator
       procedure :: get_n_par => circe1_data_get_n_par
       procedure :: get_pdg_out => circe1_data_get_pdg_out
       procedure :: get_pdg_int => circe1_data_get_pdg_int
       procedure :: allocate_sf_int => circe1_data_allocate_sf_int
       procedure :: get_beam_file => circe1_data_get_beam_file
  end type circe1_data_t

  type, extends (circe1_rng_t) :: rng_obj_t
     class(rng_t), allocatable :: rng
   contains
     procedure :: generate => rng_obj_generate
  end type rng_obj_t

  type, extends (sf_int_t) :: circe1_t
     type(circe1_data_t), pointer :: data => null ()
     real(default), dimension(2) :: x = 0
     real(default), dimension(2) :: xb= 0
     real(default) :: f = 0
     logical, dimension(2) :: continuum = .true.
     logical, dimension(2) :: peak = .true.
     type(rng_obj_t) :: rng_obj
   contains
     procedure :: type_string => circe1_type_string
     procedure :: write => circe1_write
     procedure :: init => circe1_init
     procedure :: is_generator => circe1_is_generator
     procedure :: generate_free => circe1_generate_free
     procedure :: complete_kinematics => circe1_complete_kinematics
     procedure :: inverse_kinematics => circe1_inverse_kinematics
     procedure :: apply => circe1_apply
  end type circe1_t


  interface
    module subroutine circe1_data_init &
         (data, model, pdg_in, sqrts, eps, out_photon, &
          ver, rev, acc, chat, with_radiation)
      class(circe1_data_t), intent(out) :: data
      class(model_data_t), intent(in), target :: model
      type(pdg_array_t), dimension(2), intent(in) :: pdg_in
      real(default), intent(in) :: sqrts
      real(default), intent(in) :: eps
      logical, dimension(2), intent(in) :: out_photon
      character(*), intent(in) :: acc
      integer, intent(in) :: ver, rev, chat
      logical, intent(in) :: with_radiation
    end subroutine circe1_data_init
    module subroutine circe1_data_set_generator_mode (data, rng_factory)
      class(circe1_data_t), intent(inout) :: data
      class(rng_factory_t), intent(inout), allocatable :: rng_factory
    end subroutine circe1_data_set_generator_mode
    module subroutine circe1_data_check (data)
      class(circe1_data_t), intent(in) :: data
    end subroutine circe1_data_check
    module subroutine circe1_data_write (data, unit, verbose)
      class(circe1_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine circe1_data_write
    module function circe1_data_is_generator (data) result (flag)
      class(circe1_data_t), intent(in) :: data
      logical :: flag
    end function circe1_data_is_generator
    module function circe1_data_get_n_par (data) result (n)
      class(circe1_data_t), intent(in) :: data
      integer :: n
    end function circe1_data_get_n_par
    module subroutine circe1_data_get_pdg_out (data, pdg_out)
      class(circe1_data_t), intent(in) :: data
      type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    end subroutine circe1_data_get_pdg_out
    module function circe1_data_get_pdg_int (data) result (pdg)
      class(circe1_data_t), intent(in) :: data
      integer, dimension(2) :: pdg
    end function circe1_data_get_pdg_int
    module function circe1_data_get_beam_file (data) result (file)
      class(circe1_data_t), intent(in) :: data
      type(string_t) :: file
    end function circe1_data_get_beam_file
    module subroutine rng_obj_generate (rng_obj, u)
      class(rng_obj_t), intent(inout) :: rng_obj
      real(double), intent(out) :: u
    end subroutine rng_obj_generate
    module function circe1_type_string (object) result (string)
      class(circe1_t), intent(in) :: object
      type(string_t) :: string
    end function circe1_type_string
    module subroutine circe1_write (object, unit, testflag)
      class(circe1_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine circe1_write
    module subroutine circe1_init (sf_int, data)
      class(circe1_t), intent(out) :: sf_int
      class(sf_data_t), intent(in), target :: data
    end subroutine circe1_init
    module function circe1_is_generator (sf_int) result (flag)
      class(circe1_t), intent(in) :: sf_int
      logical :: flag
    end function circe1_is_generator
    module subroutine circe1_generate_free (sf_int, r, rb,  x_free)
      class(circe1_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(inout) :: x_free
    end subroutine circe1_generate_free
    module subroutine circe1_complete_kinematics &
         (sf_int, x, xb, f, r, rb, map)
      class(circe1_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), intent(in) :: rb
      logical, intent(in) :: map
    end subroutine circe1_complete_kinematics
    module subroutine circe1_inverse_kinematics &
         (sf_int, x, xb, f, r, rb, map, set_momenta)
      class(circe1_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: r
      real(default), dimension(:), intent(out) :: rb
      logical, intent(in) :: map
      logical, intent(in), optional :: set_momenta
    end subroutine circe1_inverse_kinematics
    module subroutine circe1_apply (sf_int, scale, negative_sf, rescale, i_sub)
      class(circe1_t), intent(inout) :: sf_int
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(in), optional :: rescale
      integer, intent(in), optional :: i_sub
    end subroutine circe1_apply
  end interface

contains

  subroutine circe1_data_allocate_sf_int (data, sf_int)
    class(circe1_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (circe1_t :: sf_int)
  end subroutine circe1_data_allocate_sf_int


end module sf_circe1
