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

module sf_escan

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use pdg_arrays
  use model_data
  use flavors
  use quantum_numbers
  use state_matrices
  use polarizations
  use sf_base

  implicit none
  private

  public :: escan_data_t

  type, extends(sf_data_t) :: escan_data_t
     private
     type(flavor_t), dimension(:,:), allocatable :: flv_in
     integer, dimension(2) :: n_flv = 0
     real(default) :: norm = 1
   contains
     procedure :: init => escan_data_init
     procedure :: write => escan_data_write
     procedure :: get_n_par => escan_data_get_n_par
     procedure :: get_pdg_out => escan_data_get_pdg_out
     procedure :: allocate_sf_int => escan_data_allocate_sf_int
  end type escan_data_t

  type, extends (sf_int_t) :: escan_t
     type(escan_data_t), pointer :: data => null ()
   contains
     procedure :: type_string => escan_type_string
     procedure :: write => escan_write
     procedure :: init => escan_init
     procedure :: complete_kinematics => escan_complete_kinematics
     procedure :: recover_x => escan_recover_x
     procedure :: inverse_kinematics => escan_inverse_kinematics
     procedure :: apply => escan_apply
  end type escan_t


  interface
    module subroutine escan_data_init (data, model, pdg_in, norm)
      class(escan_data_t), intent(out) :: data
      class(model_data_t), intent(in), target :: model
      type(pdg_array_t), dimension(2), intent(in) :: pdg_in
      real(default), intent(in), optional :: norm
    end subroutine escan_data_init
    module subroutine escan_data_write (data, unit, verbose)
      class(escan_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine escan_data_write
    module function escan_data_get_n_par (data) result (n)
      class(escan_data_t), intent(in) :: data
      integer :: n
    end function escan_data_get_n_par
    module subroutine escan_data_get_pdg_out (data, pdg_out)
      class(escan_data_t), intent(in) :: data
      type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    end subroutine escan_data_get_pdg_out
    module function escan_type_string (object) result (string)
      class(escan_t), intent(in) :: object
      type(string_t) :: string
    end function escan_type_string
    module subroutine escan_write (object, unit, testflag)
      class(escan_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine escan_write
    module subroutine escan_init (sf_int, data)
      class(escan_t), intent(out) :: sf_int
      class(sf_data_t), intent(in), target :: data
    end subroutine escan_init
    module subroutine escan_complete_kinematics (sf_int, x, xb, f, r, rb, map)
      class(escan_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(out) :: f
      real(default) :: sqrt_x
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), intent(in) :: rb
      logical, intent(in) :: map
    end subroutine escan_complete_kinematics
    module subroutine escan_recover_x (sf_int, x, xb, x_free)
      class(escan_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(inout), optional :: x_free
    end subroutine escan_recover_x
    module subroutine escan_inverse_kinematics &
         (sf_int, x, xb, f, r, rb, map, set_momenta)
      class(escan_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: r
      real(default), dimension(:), intent(out) :: rb
      logical, intent(in) :: map
      logical, intent(in), optional :: set_momenta
    end subroutine escan_inverse_kinematics
    module subroutine escan_apply (sf_int, scale, negative_sf, rescale, i_sub)
      class(escan_t), intent(inout) :: sf_int
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(in), optional :: rescale
      integer, intent(in), optional :: i_sub
    end subroutine escan_apply
  end interface

contains

  subroutine escan_data_allocate_sf_int (data, sf_int)
    class(escan_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (escan_t :: sf_int)
  end subroutine escan_data_allocate_sf_int


end module sf_escan
