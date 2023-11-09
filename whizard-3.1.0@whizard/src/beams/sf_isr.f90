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

module sf_isr

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use pdg_arrays
  use model_data
  use flavors
  use sf_aux
  use sf_mappings
  use sf_base
  use electron_pdfs

  implicit none
  private

  public :: isr_data_t
  public :: isr_t

  integer, parameter :: NONE = 0
  integer, parameter :: ZERO_MASS = 1
  integer, parameter :: Q_MAX_TOO_SMALL = 2
  integer, parameter :: EPS_TOO_LARGE = 3
  integer, parameter :: INVALID_ORDER = 4
  integer, parameter :: CHARGE_MIX = 5
  integer, parameter :: CHARGE_ZERO = 6
  integer, parameter :: MASS_MIX = 7

  type, extends (sf_data_t) :: isr_data_t
     private
     class(model_data_t), pointer :: model => null ()
     type(flavor_t), dimension(:), allocatable :: flv_in
     type(qed_pdf_t) :: pdf
     real(default) :: alpha = 0
     real(default) :: q_max = 0
     real(default) :: real_mass = 0
     real(default) :: mass = 0
     real(default) :: eps = 0
     real(default) :: log = 0
     logical :: recoil = .false.
     logical :: keep_energy = .true.
     integer :: order = 3
     integer :: error = NONE
   contains
     procedure :: init => isr_data_init
     procedure :: set_order => isr_data_set_order
     procedure :: check => isr_data_check
     procedure :: write => isr_data_write
     procedure :: get_n_par => isr_data_get_n_par
     procedure :: get_pdg_out => isr_data_get_pdg_out
     procedure :: get_eps => isr_data_get_eps
     procedure :: allocate_sf_int => isr_data_allocate_sf_int
  end type isr_data_t

  type, extends (sf_int_t) :: isr_t
     private
     type(isr_data_t), pointer :: data => null ()
     real(default) :: x = 0
     real(default) :: xb= 0
   contains
     procedure :: type_string => isr_type_string
     procedure :: write => isr_write
     procedure :: set_order => isr_set_order
     procedure :: complete_kinematics => isr_complete_kinematics
     procedure :: recover_x => sf_isr_recover_x
     procedure :: inverse_kinematics => isr_inverse_kinematics
     procedure :: init => isr_init
     procedure :: apply => isr_apply
  end type isr_t


  interface
    module subroutine isr_data_init (data, model, pdg_in, alpha, q_max, &
         mass, order, recoil, keep_energy)
      class(isr_data_t), intent(out) :: data
      class(model_data_t), intent(in), target :: model
      type(pdg_array_t), intent(in) :: pdg_in
      real(default), intent(in) :: alpha
      real(default), intent(in) :: q_max
      real(default), intent(in), optional :: mass
      integer, intent(in), optional :: order
      logical, intent(in), optional :: recoil
      logical, intent(in), optional :: keep_energy
    end subroutine isr_data_init
    elemental module subroutine isr_data_set_order (data, order)
      class(isr_data_t), intent(inout) :: data
      integer, intent(in) :: order
    end subroutine isr_data_set_order
    module subroutine isr_data_check (data)
      class(isr_data_t), intent(in) :: data
    end subroutine isr_data_check
    module subroutine isr_data_write (data, unit, verbose)
      class(isr_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine isr_data_write
    module function isr_data_get_n_par (data) result (n)
      class(isr_data_t), intent(in) :: data
      integer :: n
    end function isr_data_get_n_par
    module subroutine isr_data_get_pdg_out (data, pdg_out)
      class(isr_data_t), intent(in) :: data
      type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    end subroutine isr_data_get_pdg_out
    module function isr_data_get_eps (data) result (eps)
      class(isr_data_t), intent(in) :: data
      real(default) :: eps
    end function isr_data_get_eps
    module function isr_type_string (object) result (string)
      class(isr_t), intent(in) :: object
      type(string_t) :: string
    end function isr_type_string
    module subroutine isr_write (object, unit, testflag)
      class(isr_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine isr_write
    module subroutine isr_set_order (object, order)
      class(isr_t), intent(inout) :: object
      integer, intent(in) :: order
    end subroutine isr_set_order
    module subroutine isr_complete_kinematics (sf_int, x, xb, f, r, rb, map)
      class(isr_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), intent(in) :: rb
      logical, intent(in) :: map
    end subroutine isr_complete_kinematics
    module subroutine sf_isr_recover_x (sf_int, x, xb, x_free)
      class(isr_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_isr_recover_x
    module subroutine isr_inverse_kinematics &
         (sf_int, x, xb, f, r, rb, map, set_momenta)
      class(isr_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: r
      real(default), dimension(:), intent(out) :: rb
      logical, intent(in) :: map
      logical, intent(in), optional :: set_momenta
    end subroutine isr_inverse_kinematics
    module subroutine isr_init (sf_int, data)
      class(isr_t), intent(out) :: sf_int
      class(sf_data_t), intent(in), target :: data
    end subroutine isr_init
    module subroutine isr_apply (sf_int, scale, negative_sf, rescale, i_sub)
      class(isr_t), intent(inout) :: sf_int
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(in), optional :: rescale
      integer, intent(in), optional :: i_sub
    end subroutine isr_apply
  end interface

contains

  subroutine isr_data_allocate_sf_int (data, sf_int)
    class(isr_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (isr_t :: sf_int)
  end subroutine isr_data_allocate_sf_int


end module sf_isr
