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

module sf_ewa

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use pdg_arrays
  use model_data
  use flavors
  use quantum_numbers
  use state_matrices
  use polarizations
  use interactions
  use sf_aux
  use sf_base

  implicit none
  private

  public :: ewa_data_t

  integer, parameter :: NONE = 0
  integer, parameter :: ZERO_QMIN = 1
  integer, parameter :: Q_MAX_TOO_SMALL = 2
  integer, parameter :: ZERO_XMIN = 3
  integer, parameter :: MASS_MIX = 4
  integer, parameter :: ZERO_SW = 5
  integer, parameter :: ISOSPIN_MIX = 6
  integer, parameter :: WRONG_PRT = 7
  integer, parameter :: MASS_MIX_OUT = 8
  integer, parameter :: NO_EWA = 9

  type, extends(sf_data_t) :: ewa_data_t
     private
     class(model_data_t), pointer :: model => null ()
     type(flavor_t), dimension(:), allocatable :: flv_in
     type(flavor_t), dimension(:), allocatable :: flv_out
     real(default) :: pt_max
     real(default) :: sqrts
     real(default) :: x_min
     real(default) :: x_max
     real(default) :: mass
     real(default) :: m_out
     real(default) :: q_min
     real(default) :: cv
     real(default) :: ca
     real(default) :: costhw
     real(default) :: sinthw
     real(default) :: mW
     real(default) :: mZ
     real(default) :: coeff
     logical :: mass_set = .false.
     logical :: recoil = .false.
     logical :: keep_energy = .false.
     integer :: id = 0
     integer :: error = NONE
   contains
     procedure :: init => ewa_data_init
     procedure :: set_id => ewa_set_id
     procedure :: check => ewa_data_check
     procedure :: write => ewa_data_write
     procedure :: get_n_par => ewa_data_get_n_par
     procedure :: get_pdg_out => ewa_data_get_pdg_out
     procedure :: allocate_sf_int => ewa_data_allocate_sf_int
  end type ewa_data_t

  type, extends (sf_int_t) :: ewa_t
     type(ewa_data_t), pointer :: data => null ()
     real(default) :: x  = 0
     real(default) :: xb = 0
     integer :: n_me = 0
     real(default), dimension(:), allocatable :: cv
     real(default), dimension(:), allocatable :: ca
   contains
     procedure :: type_string => ewa_type_string
     procedure :: write => ewa_write
     procedure :: init => ewa_init
     procedure :: setup_constants => ewa_setup_constants
     procedure :: complete_kinematics => ewa_complete_kinematics
     procedure :: recover_x => sf_ewa_recover_x
     procedure :: inverse_kinematics => ewa_inverse_kinematics
     procedure :: apply => ewa_apply
  end type ewa_t


  interface
    module subroutine ewa_data_init (data, model, pdg_in, x_min, pt_max, &
          sqrts, recoil, keep_energy, mass)
      class(ewa_data_t), intent(inout) :: data
      class(model_data_t), intent(in), target :: model
      type(pdg_array_t), intent(in) :: pdg_in
      real(default), intent(in) :: x_min, pt_max, sqrts
      logical, intent(in) :: recoil, keep_energy
      real(default), intent(in), optional :: mass
    end subroutine ewa_data_init
    module subroutine ewa_set_id (data, id)
      class(ewa_data_t), intent(inout) :: data
      integer, intent(in) :: id
    end subroutine ewa_set_id
    module subroutine ewa_data_check (data)
      class(ewa_data_t), intent(in) :: data
    end subroutine ewa_data_check
    module subroutine ewa_data_write (data, unit, verbose)
      class(ewa_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine ewa_data_write
    module function ewa_data_get_n_par (data) result (n)
      class(ewa_data_t), intent(in) :: data
      integer :: n
    end function ewa_data_get_n_par
    module subroutine ewa_data_get_pdg_out (data, pdg_out)
      class(ewa_data_t), intent(in) :: data
      type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    end subroutine ewa_data_get_pdg_out
    module function ewa_type_string (object) result (string)
      class(ewa_t), intent(in) :: object
      type(string_t) :: string
    end function ewa_type_string
    module subroutine ewa_write (object, unit, testflag)
      class(ewa_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine ewa_write
    module subroutine ewa_init (sf_int, data)
      class(ewa_t), intent(out) :: sf_int
      class(sf_data_t), intent(in), target :: data
    end subroutine ewa_init
    module subroutine ewa_setup_constants (sf_int)
      class(ewa_t), intent(inout), target :: sf_int
    end subroutine ewa_setup_constants
    module subroutine ewa_complete_kinematics (sf_int, x, xb, f, r, rb, map)
      class(ewa_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), intent(in) :: rb
      logical, intent(in) :: map
    end subroutine ewa_complete_kinematics
    module subroutine sf_ewa_recover_x (sf_int, x, xb, x_free)
      class(ewa_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ewa_recover_x
    module subroutine ewa_inverse_kinematics &
         (sf_int, x, xb, f, r, rb, map, set_momenta)
      class(ewa_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: r
      real(default), dimension(:), intent(out) :: rb
      logical, intent(in) :: map
      logical, intent(in), optional :: set_momenta
    end subroutine ewa_inverse_kinematics
    module subroutine ewa_apply (sf_int, scale, negative_sf, rescale, i_sub)
      class(ewa_t), intent(inout) :: sf_int
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(in), optional :: rescale
      integer, intent(in), optional :: i_sub
    end subroutine ewa_apply
  end interface

contains

  subroutine ewa_data_allocate_sf_int (data, sf_int)
    class(ewa_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (ewa_t :: sf_int)
  end subroutine ewa_data_allocate_sf_int


end module sf_ewa
