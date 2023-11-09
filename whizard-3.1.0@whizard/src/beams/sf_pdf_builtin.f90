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

module sf_pdf_builtin

  use kinds, only: default
  use kinds, only: double
  use iso_varying_string, string_t => varying_string
  use sm_qcd
  use pdg_arrays
  use model_data
  use flavors
  use polarizations
  use sf_base

  implicit none
  private

  public :: pdf_builtin_data_t
  public :: pdf_builtin_t
  public :: alpha_qcd_pdf_builtin_t

  type, extends (sf_data_t) :: pdf_builtin_data_t
     private
     integer :: id = -1
     type (string_t) :: name
     class(model_data_t), pointer :: model => null ()
     type(flavor_t) :: flv_in
     logical :: invert
     logical :: has_photon
     logical :: photon
     logical, dimension(-6:6) :: mask
     logical :: mask_photon
     logical :: hoppet_b_matching = .false.
   contains
     procedure :: init => pdf_builtin_data_init
     procedure :: set_mask => pdf_builtin_data_set_mask
     procedure :: write => pdf_builtin_data_write
     procedure :: get_n_par => pdf_builtin_data_get_n_par
     procedure :: get_pdg_out => pdf_builtin_data_get_pdg_out
     procedure :: allocate_sf_int => pdf_builtin_data_allocate_sf_int
     procedure :: get_pdf_set => pdf_builtin_data_get_pdf_set
  end type pdf_builtin_data_t

  type, extends (sf_int_t) :: pdf_builtin_t
     type(pdf_builtin_data_t), pointer :: data => null ()
     real(default) :: x = 0
     real(default) :: q = 0
  contains
    procedure :: type_string => pdf_builtin_type_string
    procedure :: write => pdf_builtin_write
    procedure :: init => pdf_builtin_init
    procedure :: complete_kinematics => pdf_builtin_complete_kinematics
    procedure :: recover_x => pdf_builtin_recover_x
    procedure :: inverse_kinematics => pdf_builtin_inverse_kinematics
    procedure :: apply => pdf_builtin_apply
  end type pdf_builtin_t

  type, extends (alpha_qcd_t) :: alpha_qcd_pdf_builtin_t
     type(string_t) :: pdfset_name
     integer :: pdfset_id = -1
   contains
     procedure :: write => alpha_qcd_pdf_builtin_write
     procedure :: get => alpha_qcd_pdf_builtin_get
     procedure :: init => alpha_qcd_pdf_builtin_init
  end type alpha_qcd_pdf_builtin_t


  interface
    module subroutine pdf_builtin_data_init (data, &
         model, pdg_in, name, path, hoppet_b_matching)
      class(pdf_builtin_data_t), intent(out) :: data
      class(model_data_t), intent(in), target :: model
      type(pdg_array_t), intent(in) :: pdg_in
      type(string_t), intent(in) :: name
      type(string_t), intent(in) :: path
      logical, intent(in), optional :: hoppet_b_matching
    end subroutine pdf_builtin_data_init
    module subroutine pdf_builtin_data_set_mask (data, mask)
      class(pdf_builtin_data_t), intent(inout) :: data
      logical, dimension(-6:6), intent(in) :: mask
    end subroutine pdf_builtin_data_set_mask
    module subroutine pdf_builtin_data_write (data, unit, verbose)
      class(pdf_builtin_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine pdf_builtin_data_write
    module function pdf_builtin_data_get_n_par (data) result (n)
      class(pdf_builtin_data_t), intent(in) :: data
      integer :: n
    end function pdf_builtin_data_get_n_par
    module subroutine pdf_builtin_data_get_pdg_out (data, pdg_out)
      class(pdf_builtin_data_t), intent(in) :: data
      type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    end subroutine pdf_builtin_data_get_pdg_out
    elemental module function pdf_builtin_data_get_pdf_set &
         (data) result (pdf_set)
      class(pdf_builtin_data_t), intent(in) :: data
      integer :: pdf_set
    end function pdf_builtin_data_get_pdf_set
    module function pdf_builtin_type_string (object) result (string)
      class(pdf_builtin_t), intent(in) :: object
      type(string_t) :: string
    end function pdf_builtin_type_string
    module subroutine pdf_builtin_write (object, unit, testflag)
      class(pdf_builtin_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine pdf_builtin_write
    module subroutine pdf_builtin_init (sf_int, data)
      class(pdf_builtin_t), intent(out) :: sf_int
      class(sf_data_t), intent(in), target :: data
    end subroutine pdf_builtin_init
    module subroutine pdf_builtin_complete_kinematics &
         (sf_int, x, xb, f, r, rb, map)
      class(pdf_builtin_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), intent(in) :: rb
      logical, intent(in) :: map
    end subroutine pdf_builtin_complete_kinematics
    module subroutine pdf_builtin_recover_x (sf_int, x, xb, x_free)
      class(pdf_builtin_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(inout), optional :: x_free
    end subroutine pdf_builtin_recover_x
    module subroutine pdf_builtin_inverse_kinematics &
         (sf_int, x, xb, f, r, rb, map, set_momenta)
      class(pdf_builtin_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: r
      real(default), dimension(:), intent(out) :: rb
      logical, intent(in) :: map
      logical, intent(in), optional :: set_momenta
    end subroutine pdf_builtin_inverse_kinematics
    module subroutine pdf_builtin_apply &
         (sf_int, scale, negative_sf, rescale, i_sub)
      class(pdf_builtin_t), intent(inout) :: sf_int
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(in), optional :: rescale
      integer, intent(in), optional :: i_sub
    end subroutine pdf_builtin_apply
    module subroutine alpha_qcd_pdf_builtin_write (object, unit)
      class(alpha_qcd_pdf_builtin_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine alpha_qcd_pdf_builtin_write
    module function alpha_qcd_pdf_builtin_get (alpha_qcd, scale) result (alpha)
      class(alpha_qcd_pdf_builtin_t), intent(in) :: alpha_qcd
      real(default), intent(in) :: scale
      real(default) :: alpha
    end function alpha_qcd_pdf_builtin_get
    module subroutine alpha_qcd_pdf_builtin_init (alpha_qcd, name, path)
      class(alpha_qcd_pdf_builtin_t), intent(out) :: alpha_qcd
      type(string_t), intent(in) :: name
      type(string_t), intent(in) :: path
    end subroutine alpha_qcd_pdf_builtin_init
  end interface

contains

  subroutine pdf_builtin_data_allocate_sf_int (data, sf_int)
    class(pdf_builtin_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (pdf_builtin_t :: sf_int)
  end subroutine pdf_builtin_data_allocate_sf_int


end module sf_pdf_builtin
