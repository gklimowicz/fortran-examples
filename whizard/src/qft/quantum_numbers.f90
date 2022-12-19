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

module quantum_numbers

  use model_data
  use helicities
  use colors
  use flavors

  implicit none
  private

  public :: quantum_numbers_t
  public :: quantum_numbers_write
  public :: quantum_numbers_get_flavor
  public :: quantum_numbers_get_color
  public :: quantum_numbers_get_helicity
  public :: quantum_numbers_eq_wo_sub
  public :: assignment(=)
  public :: quantum_numbers_are_compatible
  public :: quantum_numbers_are_physical
  public :: quantum_numbers_canonicalize_color
  public :: make_color_map
  public :: quantum_numbers_translate_color
  public :: quantum_numbers_get_max_color_value
  public :: quantum_number_array_make_color_contractions
  public :: operator(.merge.)
  public :: quantum_numbers_mask_t
  public :: quantum_numbers_mask
  public :: quantum_numbers_mask_write
  public :: any
  public :: quantum_numbers_undefined

  type :: quantum_numbers_t
     private
     type(flavor_t) :: f
     type(color_t) :: c
     type(helicity_t) :: h
     integer :: sub = 0
   contains
     generic :: init => &
        quantum_numbers_init_f, &
        quantum_numbers_init_c, &
        quantum_numbers_init_h, &
        quantum_numbers_init_fc, &
        quantum_numbers_init_fh, &
        quantum_numbers_init_ch, &
        quantum_numbers_init_fch, &
        quantum_numbers_init_fs, &
        quantum_numbers_init_fhs, &
        quantum_numbers_init_fcs, &
        quantum_numbers_init_fhcs
     procedure, private :: quantum_numbers_init_f
     procedure, private :: quantum_numbers_init_c
     procedure, private :: quantum_numbers_init_h
     procedure, private :: quantum_numbers_init_fc
     procedure, private :: quantum_numbers_init_fh
     procedure, private :: quantum_numbers_init_ch
     procedure, private :: quantum_numbers_init_fch
     procedure, private :: quantum_numbers_init_fs
     procedure, private :: quantum_numbers_init_fhs
     procedure, private :: quantum_numbers_init_fcs
     procedure, private :: quantum_numbers_init_fhcs
     procedure :: write => quantum_numbers_write_single
     procedure :: write_raw => quantum_numbers_write_raw
     procedure :: read_raw => quantum_numbers_read_raw
     procedure :: get_flavor => quantum_numbers_get_flavor
     procedure :: get_color => quantum_numbers_get_color
     procedure :: get_helicity => quantum_numbers_get_helicity
     procedure :: get_sub => quantum_numbers_get_sub
     procedure :: set_color_ghost => quantum_numbers_set_color_ghost
     procedure :: set_model => quantum_numbers_set_model
     procedure :: tag_radiated => quantum_numbers_tag_radiated
     procedure :: tag_hard_process => quantum_numbers_tag_hard_process
     procedure :: set_subtraction_index => &
          quantum_numbers_set_subtraction_index
     procedure :: get_subtraction_index => &
          quantum_numbers_get_subtraction_index
     procedure :: get_color_type => quantum_numbers_get_color_type
     procedure :: are_valid => quantum_numbers_are_valid
     procedure :: are_associated => quantum_numbers_are_associated
     procedure :: are_diagonal => quantum_numbers_are_diagonal
     procedure :: is_color_ghost => quantum_numbers_is_color_ghost
     procedure :: are_hard_process => quantum_numbers_are_hard_process
     generic :: operator(.match.) => quantum_numbers_match
     generic :: operator(.fmatch.) => quantum_numbers_match_f
     generic :: operator(.hmatch.) => quantum_numbers_match_h
     generic :: operator(.fhmatch.) => quantum_numbers_match_fh
     generic :: operator(.dhmatch.) => quantum_numbers_match_hel_diag
     generic :: operator(==) => quantum_numbers_eq
     generic :: operator(/=) => quantum_numbers_neq
     procedure, private :: quantum_numbers_match
     procedure, private :: quantum_numbers_match_f
     procedure, private :: quantum_numbers_match_h
     procedure, private :: quantum_numbers_match_fh
     procedure, private :: quantum_numbers_match_hel_diag
     procedure, private :: quantum_numbers_eq
     procedure, private :: quantum_numbers_neq
     procedure :: add_color_offset => quantum_numbers_add_color_offset
     procedure :: invert_color => quantum_numbers_invert_color
     procedure :: flip_helicity => quantum_numbers_flip_helicity
     procedure :: undefine => quantum_numbers_undefine
     procedure :: undefined => quantum_numbers_undefined0
     procedure :: are_redundant => quantum_numbers_are_redundant
  end type quantum_numbers_t

  type :: quantum_numbers_mask_t
     private
     logical :: f = .false.
     logical :: c = .false.
     logical :: cg = .false.
     logical :: h = .false.
     logical :: hd = .false.
     integer :: sub = 0
   contains
     procedure :: init => quantum_numbers_mask_init
     procedure :: write => quantum_numbers_mask_write_single
     procedure :: set_flavor => quantum_numbers_mask_set_flavor
     procedure :: set_color => quantum_numbers_mask_set_color
     procedure :: set_helicity => quantum_numbers_mask_set_helicity
     procedure :: set_sub => quantum_numbers_mask_set_sub
     procedure :: assign => quantum_numbers_mask_assign
     generic :: operator(.or.) => quantum_numbers_mask_or
     procedure, private :: quantum_numbers_mask_or
     generic :: operator(.eqv.) => quantum_numbers_mask_eqv
     generic :: operator(.neqv.) => quantum_numbers_mask_neqv
     procedure, private :: quantum_numbers_mask_eqv
     procedure, private :: quantum_numbers_mask_neqv
     procedure :: diagonal_helicity => quantum_numbers_mask_diagonal_helicity
  end type quantum_numbers_mask_t


  interface quantum_numbers_write
     module procedure quantum_numbers_write_single
     module procedure quantum_numbers_write_array
  end interface
  interface assignment(=)
     module procedure quantum_numbers_assign
  end interface

  interface make_color_map
     module procedure quantum_numbers_make_color_map
  end interface make_color_map

  interface quantum_numbers_translate_color
     module procedure quantum_numbers_translate_color0
     module procedure quantum_numbers_translate_color1
  end interface

  interface quantum_numbers_get_max_color_value
     module procedure quantum_numbers_get_max_color_value0
     module procedure quantum_numbers_get_max_color_value1
     module procedure quantum_numbers_get_max_color_value2
  end interface

  interface operator(.merge.)
     module procedure merge_quantum_numbers0
     module procedure merge_quantum_numbers1
  end interface

  interface quantum_numbers_mask_write
     module procedure quantum_numbers_mask_write_single
     module procedure quantum_numbers_mask_write_array
  end interface
  interface any
     module procedure quantum_numbers_mask_any
  end interface
  interface quantum_numbers_undefined
     module procedure quantum_numbers_undefined0
     module procedure quantum_numbers_undefined1
     module procedure quantum_numbers_undefined11
  end interface


  interface
    impure elemental module subroutine quantum_numbers_init_f (qn, flv)
      class(quantum_numbers_t), intent(out) :: qn
      type(flavor_t), intent(in) :: flv
    end subroutine quantum_numbers_init_f
    impure elemental module subroutine quantum_numbers_init_c (qn, col)
      class(quantum_numbers_t), intent(out) :: qn
      type(color_t), intent(in) :: col
    end subroutine quantum_numbers_init_c
    impure elemental module subroutine quantum_numbers_init_h (qn, hel)
      class(quantum_numbers_t), intent(out) :: qn
      type(helicity_t), intent(in) :: hel
    end subroutine quantum_numbers_init_h
    impure elemental module subroutine quantum_numbers_init_fc (qn, flv, col)
      class(quantum_numbers_t), intent(out) :: qn
      type(flavor_t), intent(in) :: flv
      type(color_t), intent(in) :: col
    end subroutine quantum_numbers_init_fc
    impure elemental module subroutine quantum_numbers_init_fh (qn, flv, hel)
      class(quantum_numbers_t), intent(out) :: qn
      type(flavor_t), intent(in) :: flv
      type(helicity_t), intent(in) :: hel
    end subroutine quantum_numbers_init_fh
    impure elemental module subroutine quantum_numbers_init_ch (qn, col, hel)
      class(quantum_numbers_t), intent(out) :: qn
      type(color_t), intent(in) :: col
      type(helicity_t), intent(in) :: hel
    end subroutine quantum_numbers_init_ch
    impure elemental module subroutine quantum_numbers_init_fch (qn, flv, col, hel)
      class(quantum_numbers_t), intent(out) :: qn
      type(flavor_t), intent(in) :: flv
      type(color_t), intent(in) :: col
      type(helicity_t), intent(in) :: hel
    end subroutine quantum_numbers_init_fch
    impure elemental module subroutine quantum_numbers_init_fs (qn, flv, sub)
      class(quantum_numbers_t), intent(out) :: qn
      type(flavor_t), intent(in) :: flv
      integer, intent(in) :: sub
    end subroutine quantum_numbers_init_fs
    impure elemental module subroutine quantum_numbers_init_fhs (qn, flv, hel, sub)
      class(quantum_numbers_t), intent(out) :: qn
      type(flavor_t), intent(in) :: flv
      type(helicity_t), intent(in) :: hel
      integer, intent(in) :: sub
    end subroutine quantum_numbers_init_fhs
    impure elemental module subroutine quantum_numbers_init_fcs (qn, flv, col, sub)
      class(quantum_numbers_t), intent(out) :: qn
      type(flavor_t), intent(in) :: flv
      type(color_t), intent(in) :: col
      integer, intent(in) :: sub
    end subroutine quantum_numbers_init_fcs
    impure elemental module subroutine quantum_numbers_init_fhcs (qn, flv, hel, col, sub)
      class(quantum_numbers_t), intent(out) :: qn
      type(flavor_t), intent(in) :: flv
      type(helicity_t), intent(in) :: hel
      type(color_t), intent(in) :: col
      integer, intent(in) :: sub
    end subroutine quantum_numbers_init_fhcs
    module subroutine quantum_numbers_write_single (qn, unit, col_verbose)
      class(quantum_numbers_t), intent(in) :: qn
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: col_verbose
    end subroutine quantum_numbers_write_single
    module subroutine quantum_numbers_write_array (qn, unit, col_verbose)
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: col_verbose
    end subroutine quantum_numbers_write_array
    module subroutine quantum_numbers_write_raw (qn, u)
      class(quantum_numbers_t), intent(in) :: qn
      integer, intent(in) :: u
    end subroutine quantum_numbers_write_raw
    module subroutine quantum_numbers_read_raw (qn, u, iostat)
      class(quantum_numbers_t), intent(out) :: qn
      integer, intent(in) :: u
      integer, intent(out), optional :: iostat
    end subroutine quantum_numbers_read_raw
    impure elemental module function quantum_numbers_get_flavor (qn) result (flv)
      type(flavor_t) :: flv
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_get_flavor
    elemental module function quantum_numbers_get_color (qn) result (col)
      type(color_t) :: col
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_get_color
    elemental module function quantum_numbers_get_helicity (qn) result (hel)
      type(helicity_t) :: hel
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_get_helicity
    elemental module function quantum_numbers_get_sub (qn) result (sub)
      integer :: sub
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_get_sub
    elemental module subroutine quantum_numbers_set_color_ghost (qn, ghost)
      class(quantum_numbers_t), intent(inout) :: qn
      logical, intent(in) :: ghost
    end subroutine quantum_numbers_set_color_ghost
    impure elemental module subroutine quantum_numbers_set_model (qn, model)
      class(quantum_numbers_t), intent(inout) :: qn
      class(model_data_t), intent(in), target :: model
    end subroutine quantum_numbers_set_model
    elemental module subroutine quantum_numbers_tag_radiated (qn)
      class(quantum_numbers_t), intent(inout) :: qn
    end subroutine quantum_numbers_tag_radiated
    elemental module subroutine quantum_numbers_tag_hard_process (qn, hard)
      class(quantum_numbers_t), intent(inout) :: qn
      logical, intent(in), optional :: hard
    end subroutine quantum_numbers_tag_hard_process
    elemental module subroutine quantum_numbers_set_subtraction_index (qn, i)
      class(quantum_numbers_t), intent(inout) :: qn
      integer, intent(in) :: i
    end subroutine quantum_numbers_set_subtraction_index
    elemental module function quantum_numbers_get_subtraction_index &
         (qn) result (sub)
      integer :: sub
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_get_subtraction_index
    elemental module function quantum_numbers_get_color_type (qn) result (color_type)
      integer :: color_type
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_get_color_type
    elemental module function quantum_numbers_are_valid (qn) result (valid)
      logical :: valid
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_are_valid
    elemental module function quantum_numbers_are_associated (qn) result (flag)
      logical :: flag
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_are_associated
    elemental module function quantum_numbers_are_diagonal (qn) result (diagonal)
      logical :: diagonal
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_are_diagonal
    elemental module function quantum_numbers_is_color_ghost (qn) result (ghost)
      logical :: ghost
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_is_color_ghost
    elemental module function quantum_numbers_are_hard_process &
         (qn) result (hard_process)
      logical :: hard_process
      class(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_are_hard_process
    elemental module function quantum_numbers_match (qn1, qn2) result (match)
      logical :: match
      class(quantum_numbers_t), intent(in) :: qn1, qn2
    end function quantum_numbers_match
    elemental module function quantum_numbers_match_f (qn1, qn2) result (match)
      logical :: match
      class(quantum_numbers_t), intent(in) :: qn1, qn2
    end function quantum_numbers_match_f
    elemental module function quantum_numbers_match_h (qn1, qn2) result (match)
      logical :: match
      class(quantum_numbers_t), intent(in) :: qn1, qn2
    end function quantum_numbers_match_h
    elemental module function quantum_numbers_match_fh (qn1, qn2) result (match)
      logical :: match
      class(quantum_numbers_t), intent(in) :: qn1, qn2
    end function quantum_numbers_match_fh
    elemental module function quantum_numbers_match_hel_diag (qn1, qn2) result (match)
      logical :: match
      class(quantum_numbers_t), intent(in) :: qn1, qn2
    end function quantum_numbers_match_hel_diag
    elemental module function quantum_numbers_eq_wo_sub (qn1, qn2) result (eq)
      logical :: eq
      type(quantum_numbers_t), intent(in) :: qn1, qn2
    end function quantum_numbers_eq_wo_sub
    elemental module function quantum_numbers_eq (qn1, qn2) result (eq)
      logical :: eq
      class(quantum_numbers_t), intent(in) :: qn1, qn2
    end function quantum_numbers_eq
    elemental module function quantum_numbers_neq (qn1, qn2) result (neq)
      logical :: neq
      class(quantum_numbers_t), intent(in) :: qn1, qn2
    end function quantum_numbers_neq
    module subroutine quantum_numbers_assign (qn_out, qn_in)
      type(quantum_numbers_t), intent(out) :: qn_out
      type(quantum_numbers_t), intent(in) :: qn_in
    end subroutine quantum_numbers_assign
    elemental module function quantum_numbers_are_compatible &
         (qn1, qn2, mask) result (flag)
      logical :: flag
      type(quantum_numbers_t), intent(in) :: qn1, qn2
      type(quantum_numbers_mask_t), intent(in) :: mask
    end function quantum_numbers_are_compatible
    elemental module function quantum_numbers_are_physical (qn, mask) result (flag)
      logical :: flag
      type(quantum_numbers_t), intent(in) :: qn
      type(quantum_numbers_mask_t), intent(in) :: mask
    end function quantum_numbers_are_physical
    module subroutine quantum_numbers_canonicalize_color (qn)
      type(quantum_numbers_t), dimension(:), intent(inout) :: qn
    end subroutine quantum_numbers_canonicalize_color
    module subroutine quantum_numbers_make_color_map (map, qn1, qn2)
      integer, dimension(:,:), intent(out), allocatable :: map
      type(quantum_numbers_t), dimension(:), intent(in) :: qn1, qn2
    end subroutine quantum_numbers_make_color_map
    module subroutine quantum_numbers_translate_color0 (qn, map, offset)
      type(quantum_numbers_t), intent(inout) :: qn
      integer, dimension(:,:), intent(in) :: map
      integer, intent(in), optional :: offset
    end subroutine quantum_numbers_translate_color0
    module subroutine quantum_numbers_translate_color1 (qn, map, offset)
      type(quantum_numbers_t), dimension(:), intent(inout) :: qn
      integer, dimension(:,:), intent(in) :: map
      integer, intent(in), optional :: offset
    end subroutine quantum_numbers_translate_color1
    pure module function quantum_numbers_get_max_color_value0 (qn) result (cmax)
      integer :: cmax
      type(quantum_numbers_t), intent(in) :: qn
    end function quantum_numbers_get_max_color_value0
    pure module function quantum_numbers_get_max_color_value1 (qn) result (cmax)
      integer :: cmax
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
    end function quantum_numbers_get_max_color_value1
    pure module function quantum_numbers_get_max_color_value2 (qn) result (cmax)
      integer :: cmax
      type(quantum_numbers_t), dimension(:,:), intent(in) :: qn
    end function quantum_numbers_get_max_color_value2
    elemental module subroutine quantum_numbers_add_color_offset (qn, offset)
      class(quantum_numbers_t), intent(inout) :: qn
      integer, intent(in) :: offset
    end subroutine quantum_numbers_add_color_offset
    module subroutine quantum_number_array_make_color_contractions (qn_in, qn_out)
      type(quantum_numbers_t), dimension(:), intent(in) :: qn_in
      type(quantum_numbers_t), dimension(:,:), intent(out), allocatable :: qn_out
    end subroutine quantum_number_array_make_color_contractions
    elemental module subroutine quantum_numbers_invert_color (qn)
      class(quantum_numbers_t), intent(inout) :: qn
    end subroutine quantum_numbers_invert_color
    elemental module subroutine quantum_numbers_flip_helicity (qn)
      class(quantum_numbers_t), intent(inout) :: qn
    end subroutine quantum_numbers_flip_helicity
    module function merge_quantum_numbers0 (qn1, qn2) result (qn3)
      type(quantum_numbers_t) :: qn3
      type(quantum_numbers_t), intent(in) :: qn1, qn2
    end function merge_quantum_numbers0
    module function merge_quantum_numbers1 (qn1, qn2) result (qn3)
      type(quantum_numbers_t), dimension(:), intent(in) :: qn1, qn2
      type(quantum_numbers_t), dimension(size(qn1)) :: qn3
    end function merge_quantum_numbers1
    elemental module function quantum_numbers_mask &
         (mask_f, mask_c, mask_h, mask_cg, mask_hd) result (mask)
      type(quantum_numbers_mask_t) :: mask
      logical, intent(in) :: mask_f, mask_c, mask_h
      logical, intent(in), optional :: mask_cg
      logical, intent(in), optional :: mask_hd
    end function quantum_numbers_mask
    elemental module subroutine quantum_numbers_mask_init &
         (mask, mask_f, mask_c, mask_h, mask_cg, mask_hd)
      class(quantum_numbers_mask_t), intent(inout) :: mask
      logical, intent(in) :: mask_f, mask_c, mask_h
      logical, intent(in), optional :: mask_cg, mask_hd
    end subroutine quantum_numbers_mask_init
    module subroutine quantum_numbers_mask_write_single (mask, unit)
      class(quantum_numbers_mask_t), intent(in) :: mask
      integer, intent(in), optional :: unit
    end subroutine quantum_numbers_mask_write_single
    module subroutine quantum_numbers_mask_write_array (mask, unit)
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: mask
      integer, intent(in), optional :: unit
    end subroutine quantum_numbers_mask_write_array
    elemental module subroutine quantum_numbers_mask_set_flavor (mask, mask_f)
      class(quantum_numbers_mask_t), intent(inout) :: mask
      logical, intent(in) :: mask_f
    end subroutine quantum_numbers_mask_set_flavor
    elemental module subroutine quantum_numbers_mask_set_color (mask, mask_c, mask_cg)
      class(quantum_numbers_mask_t), intent(inout) :: mask
      logical, intent(in) :: mask_c
      logical, intent(in), optional :: mask_cg
    end subroutine quantum_numbers_mask_set_color
    elemental module subroutine quantum_numbers_mask_set_helicity (mask, mask_h, mask_hd)
      class(quantum_numbers_mask_t), intent(inout) :: mask
      logical, intent(in) :: mask_h
      logical, intent(in), optional :: mask_hd
    end subroutine quantum_numbers_mask_set_helicity
    elemental module subroutine quantum_numbers_mask_set_sub (mask, sub)
      class(quantum_numbers_mask_t), intent(inout) :: mask
      integer, intent(in) :: sub
    end subroutine quantum_numbers_mask_set_sub
    elemental module subroutine quantum_numbers_mask_assign &
         (mask, mask_in, flavor, color, helicity)
      class(quantum_numbers_mask_t), intent(inout) :: mask
      class(quantum_numbers_mask_t), intent(in) :: mask_in
      logical, intent(in), optional :: flavor, color, helicity
    end subroutine quantum_numbers_mask_assign
    module function quantum_numbers_mask_any (mask) result (match)
      logical :: match
      type(quantum_numbers_mask_t), intent(in) :: mask
    end function quantum_numbers_mask_any
    elemental module function quantum_numbers_mask_or (mask1, mask2) result (mask)
      type(quantum_numbers_mask_t) :: mask
      class(quantum_numbers_mask_t), intent(in) :: mask1, mask2
    end function quantum_numbers_mask_or
    elemental module function quantum_numbers_mask_eqv (mask1, mask2) result (eqv)
      logical :: eqv
      class(quantum_numbers_mask_t), intent(in) :: mask1, mask2
    end function quantum_numbers_mask_eqv
    elemental module function quantum_numbers_mask_neqv (mask1, mask2) result (neqv)
      logical :: neqv
      class(quantum_numbers_mask_t), intent(in) :: mask1, mask2
    end function quantum_numbers_mask_neqv
    elemental module subroutine quantum_numbers_undefine (qn, mask)
      class(quantum_numbers_t), intent(inout) :: qn
      type(quantum_numbers_mask_t), intent(in) :: mask
    end subroutine quantum_numbers_undefine
    module function quantum_numbers_undefined0 (qn, mask) result (qn_new)
      class(quantum_numbers_t), intent(in) :: qn
      type(quantum_numbers_mask_t), intent(in) :: mask
      type(quantum_numbers_t) :: qn_new
    end function quantum_numbers_undefined0
    module function quantum_numbers_undefined1 (qn, mask) result (qn_new)
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      type(quantum_numbers_mask_t), intent(in) :: mask
      type(quantum_numbers_t), dimension(size(qn)) :: qn_new
    end function quantum_numbers_undefined1
    module function quantum_numbers_undefined11 (qn, mask) result (qn_new)
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: mask
      type(quantum_numbers_t), dimension(size(qn)) :: qn_new
    end function quantum_numbers_undefined11
    elemental module function quantum_numbers_are_redundant (qn, mask) &
         result (redundant)
      logical :: redundant
      class(quantum_numbers_t), intent(in) :: qn
      type(quantum_numbers_mask_t), intent(in) :: mask
    end function quantum_numbers_are_redundant
    elemental module function quantum_numbers_mask_diagonal_helicity (mask) &
         result (flag)
      logical :: flag
      class(quantum_numbers_mask_t), intent(in) :: mask
    end function quantum_numbers_mask_diagonal_helicity
  end interface

end module quantum_numbers
