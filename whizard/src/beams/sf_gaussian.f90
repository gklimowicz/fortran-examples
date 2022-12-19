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

module sf_gaussian

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use rng_base
  use pdg_arrays
  use model_data
  use flavors
  use quantum_numbers
  use state_matrices
  use polarizations
  use sf_base

  implicit none
  private

  public :: gaussian_data_t
  public :: gaussian_t

  type, extends(sf_data_t) :: gaussian_data_t
     private
     type(flavor_t), dimension(2) :: flv_in
     real(default), dimension(2) :: spread
     class(rng_factory_t), allocatable :: rng_factory
   contains
     procedure :: init => gaussian_data_init
     procedure :: is_generator => gaussian_data_is_generator
     procedure :: get_n_par => gaussian_data_get_n_par
     procedure :: get_pdg_out => gaussian_data_get_pdg_out
     procedure :: allocate_sf_int => gaussian_data_allocate_sf_int
     procedure :: write => gaussian_data_write
  end type gaussian_data_t

  type, extends (sf_int_t) :: gaussian_t
     type(gaussian_data_t), pointer :: data => null ()
     class(rng_t), allocatable :: rng
   contains
     procedure :: type_string => gaussian_type_string
     procedure :: write => gaussian_write
     procedure :: init => gaussian_init
     procedure :: final => sf_gaussian_final
     procedure :: is_generator => gaussian_is_generator
     procedure :: generate_free => gaussian_generate_free
     procedure :: complete_kinematics => gaussian_complete_kinematics
     procedure :: inverse_kinematics => gaussian_inverse_kinematics
     procedure :: apply => gaussian_apply
  end type gaussian_t


  interface
    module subroutine gaussian_data_init &
         (data, model, pdg_in, spread, rng_factory)
      class(gaussian_data_t), intent(out) :: data
      class(model_data_t), intent(in), target :: model
      type(pdg_array_t), dimension(2), intent(in) :: pdg_in
      real(default), dimension(2), intent(in) :: spread
      class(rng_factory_t), intent(inout), allocatable :: rng_factory
    end subroutine gaussian_data_init
    module function gaussian_data_is_generator (data) result (flag)
      class(gaussian_data_t), intent(in) :: data
      logical :: flag
    end function gaussian_data_is_generator
    module function gaussian_data_get_n_par (data) result (n)
      class(gaussian_data_t), intent(in) :: data
      integer :: n
    end function gaussian_data_get_n_par
    module subroutine gaussian_data_get_pdg_out (data, pdg_out)
      class(gaussian_data_t), intent(in) :: data
      type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    end subroutine gaussian_data_get_pdg_out
    module subroutine gaussian_data_write (data, unit, verbose)
      class(gaussian_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine gaussian_data_write
    module function gaussian_type_string (object) result (string)
      class(gaussian_t), intent(in) :: object
      type(string_t) :: string
    end function gaussian_type_string
    module subroutine gaussian_write (object, unit, testflag)
      class(gaussian_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine gaussian_write
    module subroutine gaussian_init (sf_int, data)
      class(gaussian_t), intent(out) :: sf_int
      class(sf_data_t), intent(in), target :: data
    end subroutine gaussian_init
    module subroutine sf_gaussian_final (object)
      class(gaussian_t), intent(inout) :: object
    end subroutine sf_gaussian_final
    module function gaussian_is_generator (sf_int) result (flag)
      class(gaussian_t), intent(in) :: sf_int
      logical :: flag
    end function gaussian_is_generator
    module subroutine gaussian_generate_free (sf_int, r, rb, x_free)
      class(gaussian_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(inout) :: x_free
    end subroutine gaussian_generate_free
    module subroutine gaussian_complete_kinematics &
         (sf_int, x, xb, f, r, rb, map)
      class(gaussian_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), intent(in) :: rb
      logical, intent(in) :: map
    end subroutine gaussian_complete_kinematics
    module subroutine gaussian_inverse_kinematics &
         (sf_int, x, xb, f, r, rb, map, set_momenta)
      class(gaussian_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: r
      real(default), dimension(:), intent(out) :: rb
      logical, intent(in) :: map
      logical, intent(in), optional :: set_momenta
    end subroutine gaussian_inverse_kinematics
    module subroutine gaussian_apply &
         (sf_int, scale, negative_sf, rescale, i_sub)
      class(gaussian_t), intent(inout) :: sf_int
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(in), optional :: rescale
      integer, intent(in), optional :: i_sub
    end subroutine gaussian_apply
  end interface

contains

  subroutine gaussian_data_allocate_sf_int (data, sf_int)
    class(gaussian_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (gaussian_t :: sf_int)
  end subroutine gaussian_data_allocate_sf_int


end module sf_gaussian
