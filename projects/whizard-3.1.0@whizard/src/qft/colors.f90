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

module colors

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: color_t
  public :: color_init_from_array
  public :: color_write
  public :: color_get_max_value
  public :: color_canonicalize
  public :: color_array_make_contractions
  public :: make_color_map
  public :: color_translate
  public :: compute_color_factor
  public :: count_color_loops

  type :: color_t
     private
     logical :: defined = .false.
     integer, dimension(2) :: c1 = 0, c2 = 0
     logical :: ghost = .false.
   contains
     generic :: init => &
          color_init_trivial, color_init_trivial_ghost, &
          color_init_array, color_init_array_ghost, &
          color_init_arrays, color_init_arrays_ghost
     procedure, private :: color_init_trivial
     procedure, private :: color_init_trivial_ghost
     procedure, private :: color_init_array
     procedure, private :: color_init_array_ghost
     procedure, private :: color_init_arrays
     procedure, private :: color_init_arrays_ghost
     procedure :: init_col_acl => color_init_col_acl
     procedure :: set_ghost => color_set_ghost
     procedure :: undefine => color_undefine
     procedure :: write => color_write_single
     procedure :: write_raw => color_write_raw
     procedure :: read_raw => color_read_raw
     procedure :: is_defined => color_is_defined
     procedure :: is_nonzero => color_is_nonzero
     procedure :: is_diagonal => color_is_diagonal
     procedure :: is_ghost => color_is_ghost
     procedure :: get_type => color_get_type
     procedure, private :: get_number_of_indices => color_get_number_of_indices
     procedure :: get_col => color_get_col
     procedure :: get_acl => color_get_acl
     generic :: operator(.match.) => color_match
     generic :: operator(==) => color_eq
     generic :: operator(/=) => color_neq
     procedure, private ::  color_match
     procedure, private ::  color_eq
     procedure, private ::  color_neq
     procedure :: add_offset => color_add_offset
     procedure :: invert => color_invert
     generic :: operator(.merge.) => merge_colors
     procedure, private ::  merge_colors
     generic :: operator (.fuse.) => color_fusion
     procedure, private :: color_fusion
  end type color_t

  type :: entry_t
     integer, dimension(:), allocatable :: map
     type(color_t), dimension(:), allocatable :: col
     type(entry_t), pointer :: next => null ()
     logical :: nlo_event = .false.
  end type entry_t
  type :: list_t
     integer :: n = 0
     type(entry_t), pointer :: first => null ()
     type(entry_t), pointer :: last => null ()
  end type list_t


  interface color_init_from_array
     module procedure color_init_from_array1
     module procedure color_init_from_array1g
     module procedure color_init_from_array2
     module procedure color_init_from_array2g
  end interface color_init_from_array

  interface color_write
     module procedure color_write_single
     module procedure color_write_array
  end interface color_write

  interface color_get_max_value
     module procedure color_get_max_value0
     module procedure color_get_max_value1
     module procedure color_get_max_value2
  end interface color_get_max_value

  interface make_color_map
     module procedure color_make_color_map
  end interface make_color_map

  interface color_translate
     module procedure color_translate0
     module procedure color_translate0_offset
     module procedure color_translate1
  end interface color_translate


  interface
    pure module subroutine color_init_trivial (col)
      class(color_t), intent(inout) :: col
    end subroutine color_init_trivial
    pure module subroutine color_init_trivial_ghost (col, ghost)
      class(color_t), intent(inout) :: col
      logical, intent(in) :: ghost
    end subroutine color_init_trivial_ghost
    pure module subroutine color_init_array (col, c1)
      class(color_t), intent(inout) :: col
      integer, dimension(:), intent(in) :: c1
    end subroutine color_init_array
    pure module subroutine color_init_array_ghost (col, c1, ghost)
      class(color_t), intent(inout) :: col
      integer, dimension(:), intent(in) :: c1
      logical, intent(in) :: ghost
    end subroutine color_init_array_ghost
    pure module subroutine color_init_arrays (col, c1, c2)
      class(color_t), intent(inout) :: col
      integer, dimension(:), intent(in) :: c1, c2
    end subroutine color_init_arrays
    pure module subroutine color_init_arrays_ghost (col, c1, c2, ghost)
      class(color_t), intent(inout) :: col
      integer, dimension(:), intent(in) :: c1, c2
      logical, intent(in) :: ghost
    end subroutine color_init_arrays_ghost
    elemental module subroutine color_init_col_acl (col, col_in, acl_in)
      class(color_t), intent(inout) :: col
      integer, intent(in) :: col_in, acl_in
    end subroutine color_init_col_acl
    pure module subroutine color_init_from_array1 (col, c1)
      type(color_t), intent(inout) :: col
      integer, dimension(:), intent(in) :: c1
    end subroutine color_init_from_array1
    pure module subroutine color_init_from_array1g (col, c1, ghost)
      type(color_t), intent(inout) :: col
      integer, dimension(:), intent(in) :: c1
      logical, intent(in) :: ghost
    end subroutine color_init_from_array1g
    pure module subroutine color_init_from_array2 (col, c1)
      integer, dimension(:,:), intent(in) :: c1
      type(color_t), dimension(:), intent(inout) :: col
    end subroutine color_init_from_array2
    pure module subroutine color_init_from_array2g (col, c1, ghost)
      integer, dimension(:,:), intent(in) :: c1
      type(color_t), dimension(:), intent(out) :: col
      logical, intent(in), dimension(:) :: ghost
    end subroutine color_init_from_array2g
    elemental module subroutine color_set_ghost (col, ghost)
      class(color_t), intent(inout) :: col
      logical, intent(in) :: ghost
    end subroutine color_set_ghost
    elemental module subroutine color_undefine (col, undefine_ghost)
      class(color_t), intent(inout) :: col
      logical, intent(in), optional :: undefine_ghost
    end subroutine color_undefine
    module subroutine color_write_single (col, unit)
      class(color_t), intent(in) :: col
      integer, intent(in), optional :: unit
    end subroutine color_write_single
    module subroutine color_write_array (col, unit)
      type(color_t), dimension(:), intent(in) :: col
      integer, intent(in), optional :: unit
    end subroutine color_write_array
    module subroutine color_write_raw (col, u)
      class(color_t), intent(in) :: col
      integer, intent(in) :: u
    end subroutine color_write_raw
    module subroutine color_read_raw (col, u, iostat)
      class(color_t), intent(inout) :: col
      integer, intent(in) :: u
      integer, intent(out), optional :: iostat
    end subroutine color_read_raw
    elemental module function color_is_defined (col) result (defined)
      logical :: defined
      class(color_t), intent(in) :: col
    end function color_is_defined
    elemental module function color_is_nonzero (col) result (flag)
      logical :: flag
      class(color_t), intent(in) :: col
    end function color_is_nonzero
    elemental module function color_is_diagonal (col) result (diagonal)
      logical :: diagonal
      class(color_t), intent(in) :: col
    end function color_is_diagonal
    elemental module function color_is_ghost (col) result (ghost)
      logical :: ghost
      class(color_t), intent(in) :: col
    end function color_is_ghost
    elemental module function color_get_type (col) result (ctype)
      class(color_t), intent(in) :: col
      integer :: ctype
    end function color_get_type
    elemental module function color_get_number_of_indices (col) result (n)
      integer :: n
      class(color_t), intent(in) :: col
    end function color_get_number_of_indices
    elemental module function color_get_col (col) result (c)
      integer :: c
      class(color_t), intent(in) :: col
    end function color_get_col
    elemental module function color_get_acl (col) result (c)
      integer :: c
      class(color_t), intent(in) :: col
    end function color_get_acl
    elemental module function color_get_max_value0 (col) result (cmax)
      integer :: cmax
      type(color_t), intent(in) :: col
    end function color_get_max_value0
    pure module function color_get_max_value1 (col) result (cmax)
      integer :: cmax
      type(color_t), dimension(:), intent(in) :: col
    end function color_get_max_value1
    pure module function color_get_max_value2 (col) result (cmax)
      integer :: cmax
      type(color_t), dimension(:,:), intent(in) :: col
    end function color_get_max_value2
    elemental module function color_match (col1, col2) result (eq)
      logical :: eq
      class(color_t), intent(in) :: col1, col2
    end function color_match
    elemental module function color_eq (col1, col2) result (eq)
      logical :: eq
      class(color_t), intent(in) :: col1, col2
    end function color_eq
    elemental module function color_neq (col1, col2) result (neq)
      logical :: neq
      class(color_t), intent(in) :: col1, col2
    end function color_neq
    elemental module subroutine color_add_offset (col, offset)
      class(color_t), intent(inout) :: col
      integer, intent(in) :: offset
    end subroutine color_add_offset
    module subroutine color_canonicalize (col)
      type(color_t), dimension(:), intent(inout) :: col
    end subroutine color_canonicalize
    module subroutine color_array_make_contractions (col_in, col_out)
      type(color_t), dimension(:), intent(in) :: col_in
      type(color_t), dimension(:,:), intent(out), allocatable :: col_out
    end subroutine color_array_make_contractions
    elemental module subroutine color_invert (col)
      class(color_t), intent(inout) :: col
    end subroutine color_invert
    module subroutine color_make_color_map (map, col1, col2)
      integer, dimension(:,:), intent(out), allocatable :: map
      type(color_t), dimension(:), intent(in) :: col1, col2
    end subroutine color_make_color_map
    module subroutine color_translate0 (col, map)
      type(color_t), intent(inout) :: col
      integer, dimension(:,:), intent(in) :: map
    end subroutine color_translate0
    module subroutine color_translate0_offset (col, map, offset)
      type(color_t), intent(inout) :: col
      integer, dimension(:,:), intent(in) :: map
      integer, intent(in) :: offset
    end subroutine color_translate0_offset
    module subroutine color_translate1 (col, map, offset)
      type(color_t), dimension(:), intent(inout) :: col
      integer, dimension(:,:), intent(in) :: map
      integer, intent(in), optional :: offset
    end subroutine color_translate1
    elemental module function merge_colors (col1, col2) result (col)
      type(color_t) :: col
      class(color_t), intent(in) :: col1, col2
    end function merge_colors
    module function color_fusion (col1, col2) result (col)
      class(color_t), intent(in) :: col1, col2
      type(color_t) :: col
    end function color_fusion
    module function compute_color_factor (col1, col2, nc) result (factor)
      real(default) :: factor
      type(color_t), dimension(:), intent(in) :: col1, col2
      integer, intent(in), optional :: nc
    end function compute_color_factor
    module function count_color_loops (col) result (count)
      integer :: count
      type(color_t), dimension(:), intent(in) :: col
    end function count_color_loops
  end interface

end module colors
