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

module process_constants

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use pdg_arrays

  implicit none
  private

  public :: process_constants_t

  type :: process_constants_t
     type(string_t) :: id
     type(string_t) :: model_name
     character(32) :: md5sum = ""
     logical :: openmp_supported = .false.
     integer :: n_in  = 0
     integer :: n_out = 0
     integer :: n_flv = 0
     integer :: n_hel = 0
     integer :: n_col = 0
     integer :: n_cin = 0
     integer :: n_cf  = 0
     integer, dimension(:,:), allocatable :: flv_state
     integer, dimension(:,:), allocatable :: hel_state
     integer, dimension(:,:,:), allocatable :: col_state
     logical, dimension(:,:), allocatable :: ghost_flag
     complex(default), dimension(:), allocatable :: color_factors
     integer, dimension(:,:), allocatable :: cf_index
     integer, dimension(:), allocatable :: eqv_flv_index
     integer, dimension(:), allocatable :: eqv_hel_index
  contains
    procedure :: get_n_tot => process_constants_get_n_tot
    procedure :: get_flv_state => process_constants_get_flv_state
    procedure :: get_n_flv => process_constants_get_n_flv
    procedure :: get_n_hel => process_constants_get_n_hel
    procedure :: get_hel_state => process_constants_get_hel_state
    procedure :: get_col_state => process_constants_get_col_state
    procedure :: get_ghost_flag => process_constants_get_ghost_flag
    procedure :: get_color_factors => process_constants_get_color_factors
    procedure :: get_cf_index => process_constants_get_cf_index
    procedure :: set_flv_state => process_constants_set_flv_state
    procedure :: set_col_state => process_constants_set_col_state
    procedure :: set_cf_index => process_constants_set_cf_index
    procedure :: set_color_factors => process_constants_set_color_factors
    procedure :: set_ghost_flag => process_constants_set_ghost_flag
    procedure :: get_pdg_in => process_constants_get_pdg_in
    procedure :: compute_md5sum => process_constants_compute_md5sum
    procedure :: fill_unit_for_md5sum => process_constants_fill_unit_for_md5sum
    procedure :: write => process_constants_write
  end type process_constants_t


  interface
    elemental module function process_constants_get_n_tot (prc_const) result (n_tot)
      integer :: n_tot
      class(process_constants_t), intent(in) :: prc_const
    end function process_constants_get_n_tot
    module subroutine process_constants_get_flv_state (prc_const, flv_state)
      class(process_constants_t), intent(in) :: prc_const
      integer, dimension(:,:), allocatable, intent(out) :: flv_state
    end subroutine process_constants_get_flv_state
    module function process_constants_get_n_flv (data) result (n_flv)
      integer :: n_flv
      class(process_constants_t), intent(in) :: data
    end function process_constants_get_n_flv
    module function process_constants_get_n_hel (data) result (n_hel)
      integer :: n_hel
      class(process_constants_t), intent(in) :: data
    end function process_constants_get_n_hel
    module subroutine process_constants_get_hel_state (prc_const, hel_state)
      class(process_constants_t), intent(in) :: prc_const
      integer, dimension(:,:), allocatable, intent(out) :: hel_state
    end subroutine process_constants_get_hel_state
    module subroutine process_constants_get_col_state (prc_const, col_state)
      class(process_constants_t), intent(in) :: prc_const
      integer, dimension(:,:,:), allocatable, intent(out) :: col_state
    end subroutine process_constants_get_col_state
    module subroutine process_constants_get_ghost_flag (prc_const, ghost_flag)
      class(process_constants_t), intent(in) :: prc_const
      logical, dimension(:,:), allocatable, intent(out) :: ghost_flag
    end subroutine process_constants_get_ghost_flag
    module subroutine process_constants_get_color_factors (prc_const, col_facts)
      class(process_constants_t), intent(in) :: prc_const
      complex(default), dimension(:), allocatable, intent(out) :: col_facts
    end subroutine process_constants_get_color_factors
    module subroutine process_constants_get_cf_index (prc_const, cf_index)
      class(process_constants_t), intent(in) :: prc_const
      integer, intent(out), dimension(:,:), allocatable :: cf_index
    end subroutine process_constants_get_cf_index
    module subroutine process_constants_set_flv_state (prc_const, flv_state)
      class(process_constants_t), intent(inout) :: prc_const
      integer, intent(in), dimension(:,:), allocatable :: flv_state
    end subroutine process_constants_set_flv_state
    module subroutine process_constants_set_col_state (prc_const, col_state)
      class(process_constants_t), intent(inout) :: prc_const
      integer, intent(in), dimension(:,:,:), allocatable :: col_state
    end subroutine process_constants_set_col_state
    module subroutine process_constants_set_cf_index (prc_const, cf_index)
      class(process_constants_t), intent(inout) :: prc_const
      integer, dimension(:,:), intent(in), allocatable :: cf_index
    end subroutine process_constants_set_cf_index
    module subroutine process_constants_set_color_factors (prc_const, color_factors)
      class(process_constants_t), intent(inout) :: prc_const
      complex(default), dimension(:), intent(in), allocatable :: color_factors
    end subroutine process_constants_set_color_factors
    module subroutine process_constants_set_ghost_flag (prc_const, ghost_flag)
      class(process_constants_t), intent(inout) :: prc_const
      logical, dimension(:,:), allocatable, intent(in) :: ghost_flag
    end subroutine process_constants_set_ghost_flag
    module function process_constants_get_pdg_in (prc_const) result (pdg_in)
      type(pdg_array_t), dimension(:), allocatable :: pdg_in
      class(process_constants_t), intent(in) :: prc_const
    end function process_constants_get_pdg_in
    module subroutine process_constants_compute_md5sum (prc_const, include_id)
      class(process_constants_t), intent(inout) :: prc_const
      logical, intent(in) :: include_id
    end subroutine process_constants_compute_md5sum
    module function process_constants_fill_unit_for_md5sum &
         (prc_const, include_id) result (unit)
      integer :: unit
      class(process_constants_t), intent(in) :: prc_const
      logical, intent(in) :: include_id
    end function process_constants_fill_unit_for_md5sum
    module subroutine process_constants_write (prc_const, unit)
      class(process_constants_t), intent(in) :: prc_const
      integer, intent(in), optional :: unit
    end subroutine process_constants_write
  end interface

end module process_constants
