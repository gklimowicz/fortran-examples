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

module sf_epa

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
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

  public :: EPA_MODE_DEFAULT
  public :: EPA_MODE_BUDNEV_617
  public :: EPA_MODE_BUDNEV_616E
  public :: EPA_MODE_LOG_POWER
  public :: EPA_MODE_LOG_SIMPLE
  public :: EPA_MODE_LOG
  public :: epa_data_t

  integer, parameter :: EPA_MODE_DEFAULT = 0
  integer, parameter :: EPA_MODE_BUDNEV_617 = 0
  integer, parameter :: EPA_MODE_BUDNEV_616E = 1
  integer, parameter :: EPA_MODE_LOG_POWER = 2
  integer, parameter :: EPA_MODE_LOG_SIMPLE = 3
  integer, parameter :: EPA_MODE_LOG = 4

  integer, parameter :: NONE = 0
  integer, parameter :: ZERO_QMIN = 1
  integer, parameter :: Q_MAX_TOO_SMALL = 2
  integer, parameter :: ZERO_XMIN = 3
  integer, parameter :: MASS_MIX = 4
  integer, parameter :: NO_EPA = 5

  type, extends(sf_data_t) :: epa_data_t
     private
     class(model_data_t), pointer :: model => null ()
     type(flavor_t), dimension(:), allocatable :: flv_in
     real(default) :: alpha
     real(default) :: x_min
     real(default) :: x_max
     real(default) :: q_min
     real(default) :: q_max
     real(default) :: E_max
     real(default) :: mass
     real(default) :: log
     real(default) :: a
     real(default) :: c0
     real(default) :: c1
     real(default) :: dc
     integer :: mode = EPA_MODE_DEFAULT
     integer :: error = NONE
     logical :: recoil = .false.
     logical :: keep_energy = .true.
   contains
     procedure :: init => epa_data_init
     procedure :: check => epa_data_check
     procedure :: write => epa_data_write
     procedure :: get_n_par => epa_data_get_n_par
     procedure :: get_pdg_out => epa_data_get_pdg_out
     procedure :: allocate_sf_int => epa_data_allocate_sf_int
  end type epa_data_t

  type, extends (sf_int_t) :: epa_t
     type(epa_data_t), pointer :: data => null ()
     real(default) :: x  = 0
     real(default) :: xb = 0
     real(default) :: E  = 0
     real(default), dimension(:), allocatable :: charge2
   contains
     procedure :: type_string => epa_type_string
     procedure :: write => epa_write
     procedure :: init => epa_init
     procedure :: setup_constants => epa_setup_constants
     procedure :: complete_kinematics => epa_complete_kinematics
     procedure :: recover_x => sf_epa_recover_x
     procedure :: inverse_kinematics => epa_inverse_kinematics
     procedure :: apply => epa_apply
  end type epa_t


  interface
    module subroutine epa_data_init (data, model, mode, pdg_in, alpha, &
         x_min, q_min, q_max, mass, recoil, keep_energy)
      class(epa_data_t), intent(inout) :: data
      class(model_data_t), intent(in), target :: model
      type(pdg_array_t), intent(in) :: pdg_in
      integer, intent(in) :: mode
      real(default), intent(in) :: alpha, x_min, q_min, q_max
      real(default), intent(in), optional :: mass
      logical, intent(in), optional :: recoil
      logical, intent(in), optional :: keep_energy
    end subroutine epa_data_init
    module subroutine epa_data_check (data)
      class(epa_data_t), intent(in) :: data
    end subroutine epa_data_check
    module subroutine epa_data_write (data, unit, verbose)
      class(epa_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine epa_data_write
    module function epa_data_get_n_par (data) result (n)
      class(epa_data_t), intent(in) :: data
      integer :: n
    end function epa_data_get_n_par
    module subroutine epa_data_get_pdg_out (data, pdg_out)
      class(epa_data_t), intent(in) :: data
      type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    end subroutine epa_data_get_pdg_out
    module function epa_type_string (object) result (string)
      class(epa_t), intent(in) :: object
      type(string_t) :: string
    end function epa_type_string
    module subroutine epa_write (object, unit, testflag)
      class(epa_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine epa_write
    module subroutine epa_init (sf_int, data)
      class(epa_t), intent(out) :: sf_int
      class(sf_data_t), intent(in), target :: data
    end subroutine epa_init
    module subroutine epa_setup_constants (sf_int)
      class(epa_t), intent(inout), target :: sf_int
    end subroutine epa_setup_constants
    module subroutine epa_complete_kinematics (sf_int, x, xb, f, r, rb, map)
      class(epa_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), intent(in) :: rb
      logical, intent(in) :: map
    end subroutine epa_complete_kinematics
    module subroutine sf_epa_recover_x (sf_int, x, xb, x_free)
      class(epa_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_epa_recover_x
    module subroutine epa_inverse_kinematics &
         (sf_int, x, xb, f, r, rb, map, set_momenta)
      class(epa_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: r
      real(default), dimension(:), intent(out) :: rb
      logical, intent(in) :: map
      logical, intent(in), optional :: set_momenta
    end subroutine epa_inverse_kinematics
    module subroutine epa_apply (sf_int, scale, negative_sf, rescale, i_sub)
      class(epa_t), intent(inout) :: sf_int
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(in), optional :: rescale
      integer, intent(in), optional :: i_sub
    end subroutine epa_apply
  end interface

contains

  subroutine epa_data_allocate_sf_int (data, sf_int)
    class(epa_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (epa_t :: sf_int)
  end subroutine epa_data_allocate_sf_int


end module sf_epa
