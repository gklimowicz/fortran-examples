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

module sf_circe2

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use rng_base
  use selectors
  use pdg_arrays
  use model_data
  use flavors
  use polarizations
  use sf_base
  use circe2, circe2_rng_t => rng_type !NODEP!

  implicit none
  private

  public :: circe2_data_t
  public :: circe2_t

  type, extends (sf_data_t) :: circe2_data_t
     private
     class(model_data_t), pointer :: model => null ()
     type(flavor_t), dimension(2) :: flv_in
     integer, dimension(2) :: pdg_in
     real(default) :: sqrts = 0
     logical :: polarized = .false.
     logical :: beams_polarized = .false.
     class(rng_factory_t), allocatable :: rng_factory
     type(string_t) :: filename
     type(string_t) :: file
     type(string_t) :: design
     real(default) :: lumi = 0
     real(default), dimension(4) :: lumi_hel_frac = 0
     integer, dimension(0:4) :: h1 = [0, -1, -1, 1, 1]
     integer, dimension(0:4) :: h2 = [0, -1,  1,-1, 1]
     integer :: error = 1
   contains
       procedure :: init => circe2_data_init
       procedure :: set_generator_mode => circe2_data_set_generator_mode
       procedure :: check_file => circe2_check_file
       procedure :: check => circe2_data_check
       procedure :: write => circe2_data_write
       procedure :: is_generator => circe2_data_is_generator
       procedure :: get_n_par => circe2_data_get_n_par
       procedure :: get_pdg_out => circe2_data_get_pdg_out
       procedure :: allocate_sf_int => circe2_data_allocate_sf_int
       procedure :: get_beam_file => circe2_data_get_beam_file
  end type circe2_data_t

  type(circe2_state) :: circe2_global_state

  type, extends (circe2_rng_t) :: rng_obj_t
     class(rng_t), allocatable :: rng
   contains
     procedure :: generate => rng_obj_generate
  end type rng_obj_t

  type, extends (sf_int_t) :: circe2_t
     type(circe2_data_t), pointer :: data => null ()
     type(rng_obj_t) :: rng_obj
     type(selector_t) :: selector
     integer :: h_sel = 0
   contains
     procedure :: type_string => circe2_type_string
     procedure :: write => circe2_write
     procedure :: init => circe2_init
     procedure :: is_generator => circe2_is_generator
     procedure :: generate_free => circe2_generate_whizard_free
     procedure :: complete_kinematics => circe2_complete_kinematics
     procedure :: inverse_kinematics => circe2_inverse_kinematics
     procedure :: apply => circe2_apply
  end type circe2_t


  interface
    module subroutine circe2_data_init (data, os_data, model, pdg_in, &
         sqrts, polarized, beam_pol, file, design)
      class(circe2_data_t), intent(out) :: data
      type(os_data_t), intent(in) :: os_data
      class(model_data_t), intent(in), target :: model
      type(pdg_array_t), dimension(2), intent(in) :: pdg_in
      real(default), intent(in) :: sqrts
      logical, intent(in) :: polarized, beam_pol
      type(string_t), intent(in) :: file, design
    end subroutine circe2_data_init
    module subroutine circe2_data_set_generator_mode (data, rng_factory)
      class(circe2_data_t), intent(inout) :: data
      class(rng_factory_t), intent(inout), allocatable :: rng_factory
    end subroutine circe2_data_set_generator_mode
    module subroutine circe2_check_file (data, os_data)
      class(circe2_data_t), intent(inout) :: data
      type(os_data_t), intent(in) :: os_data
    end subroutine circe2_check_file
    module subroutine circe2_data_check (data)
      class(circe2_data_t), intent(in) :: data
    end subroutine circe2_data_check
    module subroutine circe2_data_write (data, unit, verbose)
      class(circe2_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine circe2_data_write
    module function circe2_data_is_generator (data) result (flag)
      class(circe2_data_t), intent(in) :: data
      logical :: flag
    end function circe2_data_is_generator
    module function circe2_data_get_n_par (data) result (n)
      class(circe2_data_t), intent(in) :: data
      integer :: n
    end function circe2_data_get_n_par
    module subroutine circe2_data_get_pdg_out (data, pdg_out)
      class(circe2_data_t), intent(in) :: data
      type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    end subroutine circe2_data_get_pdg_out
    module function circe2_data_get_beam_file (data) result (file)
      class(circe2_data_t), intent(in) :: data
      type(string_t) :: file
    end function circe2_data_get_beam_file
    module subroutine rng_obj_generate (rng_obj, u)
      class(rng_obj_t), intent(inout) :: rng_obj
      real(default), intent(out) :: u
    end subroutine rng_obj_generate
    module function circe2_type_string (object) result (string)
      class(circe2_t), intent(in) :: object
      type(string_t) :: string
    end function circe2_type_string
    module subroutine circe2_write (object, unit, testflag)
      class(circe2_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine circe2_write
    module subroutine circe2_init (sf_int, data)
      class(circe2_t), intent(out) :: sf_int
      class(sf_data_t), intent(in), target :: data
    end subroutine circe2_init
    module function circe2_is_generator (sf_int) result (flag)
      class(circe2_t), intent(in) :: sf_int
      logical :: flag
    end function circe2_is_generator
    module subroutine circe2_generate_whizard_free (sf_int, r, rb, x_free)
      class(circe2_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(inout) :: x_free
    end subroutine circe2_generate_whizard_free
    module subroutine circe2_generate_whizard (x, pdg, hel, rng_obj)
      real(default), dimension(2), intent(out) :: x
      integer, dimension(2), intent(in) :: pdg
      integer, dimension(2), intent(in) :: hel
      class(rng_obj_t), intent(inout) :: rng_obj
    end subroutine circe2_generate_whizard
    module subroutine circe2_complete_kinematics (sf_int, x, xb, f, r, rb, map)
      class(circe2_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), intent(in) :: rb
      logical, intent(in) :: map
    end subroutine circe2_complete_kinematics
    module subroutine circe2_inverse_kinematics &
         (sf_int, x, xb, f, r, rb, map, set_momenta)
      class(circe2_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: r
      real(default), dimension(:), intent(out) :: rb
      logical, intent(in) :: map
      logical, intent(in), optional :: set_momenta
    end subroutine circe2_inverse_kinematics
    module subroutine circe2_apply (sf_int, scale, negative_sf, rescale, i_sub)
      class(circe2_t), intent(inout) :: sf_int
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(in), optional :: rescale
      integer, intent(in), optional :: i_sub
    end subroutine circe2_apply
  end interface

contains

  subroutine circe2_data_allocate_sf_int (data, sf_int)
    class(circe2_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (circe2_t :: sf_int)
  end subroutine circe2_data_allocate_sf_int


end module sf_circe2
