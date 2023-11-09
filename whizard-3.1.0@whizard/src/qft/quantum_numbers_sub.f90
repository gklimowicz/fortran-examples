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

submodule (quantum_numbers) quantum_numbers_s

  use io_units

  implicit none

contains

  impure elemental module subroutine quantum_numbers_init_f (qn, flv)
    class(quantum_numbers_t), intent(out) :: qn
    type(flavor_t), intent(in) :: flv
    qn%f = flv
    call qn%c%undefine ()
    call qn%h%undefine ()
    qn%sub = 0
  end subroutine quantum_numbers_init_f

  impure elemental module subroutine quantum_numbers_init_c (qn, col)
    class(quantum_numbers_t), intent(out) :: qn
    type(color_t), intent(in) :: col
    call qn%f%undefine ()
    qn%c = col
    call qn%h%undefine ()
    qn%sub = 0
  end subroutine quantum_numbers_init_c

  impure elemental module subroutine quantum_numbers_init_h (qn, hel)
    class(quantum_numbers_t), intent(out) :: qn
    type(helicity_t), intent(in) :: hel
    call qn%f%undefine ()
    call qn%c%undefine ()
    qn%h = hel
    qn%sub = 0
  end subroutine quantum_numbers_init_h

  impure elemental module subroutine quantum_numbers_init_fc (qn, flv, col)
    class(quantum_numbers_t), intent(out) :: qn
    type(flavor_t), intent(in) :: flv
    type(color_t), intent(in) :: col
    qn%f = flv
    qn%c = col
    call qn%h%undefine ()
    qn%sub = 0
  end subroutine quantum_numbers_init_fc

  impure elemental module subroutine quantum_numbers_init_fh (qn, flv, hel)
    class(quantum_numbers_t), intent(out) :: qn
    type(flavor_t), intent(in) :: flv
    type(helicity_t), intent(in) :: hel
    qn%f = flv
    call qn%c%undefine ()
    qn%h = hel
    qn%sub = 0
  end subroutine quantum_numbers_init_fh

  impure elemental module subroutine quantum_numbers_init_ch (qn, col, hel)
    class(quantum_numbers_t), intent(out) :: qn
    type(color_t), intent(in) :: col
    type(helicity_t), intent(in) :: hel
    call qn%f%undefine ()
    qn%c = col
    qn%h = hel
    qn%sub = 0
  end subroutine quantum_numbers_init_ch

  impure elemental module subroutine quantum_numbers_init_fch (qn, flv, col, hel)
    class(quantum_numbers_t), intent(out) :: qn
    type(flavor_t), intent(in) :: flv
    type(color_t), intent(in) :: col
    type(helicity_t), intent(in) :: hel
    qn%f = flv
    qn%c = col
    qn%h = hel
    qn%sub = 0
  end subroutine quantum_numbers_init_fch

  impure elemental module subroutine quantum_numbers_init_fs (qn, flv, sub)
    class(quantum_numbers_t), intent(out) :: qn
    type(flavor_t), intent(in) :: flv
    integer, intent(in) :: sub
    qn%f = flv; qn%sub = sub
  end subroutine quantum_numbers_init_fs

  impure elemental module subroutine quantum_numbers_init_fhs (qn, flv, hel, sub)
    class(quantum_numbers_t), intent(out) :: qn
    type(flavor_t), intent(in) :: flv
    type(helicity_t), intent(in) :: hel
    integer, intent(in) :: sub
    qn%f = flv; qn%h = hel; qn%sub = sub
  end subroutine quantum_numbers_init_fhs

  impure elemental module subroutine quantum_numbers_init_fcs (qn, flv, col, sub)
    class(quantum_numbers_t), intent(out) :: qn
    type(flavor_t), intent(in) :: flv
    type(color_t), intent(in) :: col
    integer, intent(in) :: sub
    qn%f = flv; qn%c = col; qn%sub = sub
  end subroutine quantum_numbers_init_fcs

  impure elemental module subroutine quantum_numbers_init_fhcs (qn, flv, hel, col, sub)
    class(quantum_numbers_t), intent(out) :: qn
    type(flavor_t), intent(in) :: flv
    type(helicity_t), intent(in) :: hel
    type(color_t), intent(in) :: col
    integer, intent(in) :: sub
    qn%f = flv; qn%h = hel; qn%c = col; qn%sub = sub
  end subroutine quantum_numbers_init_fhcs

  module subroutine quantum_numbers_write_single (qn, unit, col_verbose)
    class(quantum_numbers_t), intent(in) :: qn
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: col_verbose
    integer :: u
    logical :: col_verb
    u = given_output_unit (unit);  if (u < 0)  return
    col_verb = .false.;  if (present (col_verbose)) col_verb = col_verbose
    write (u, "(A)", advance = "no")  "["
    if (qn%f%is_defined ()) then
       call qn%f%write (u)
       if (qn%c%is_nonzero () .or. qn%h%is_defined ()) &
            write (u, "(1x)", advance = "no")
    end if
    if (col_verb) then
       if (qn%c%is_defined () .or. qn%c%is_ghost ()) then
          call color_write (qn%c, u)
          if (qn%h%is_defined ())  write (u, "(1x)", advance = "no")
       end if
    else
       if (qn%c%is_nonzero () .or. qn%c%is_ghost ()) then
          call color_write (qn%c, u)
          if (qn%h%is_defined ())  write (u, "(1x)", advance = "no")
       end if
    end if
    if (qn%h%is_defined ()) then
       call qn%h%write (u)
    end if
    if (qn%sub > 0) &
       write (u, "(A,I0)", advance = "no") " SUB = ", qn%sub
    write (u, "(A)", advance="no")  "]"
  end subroutine quantum_numbers_write_single

  module subroutine quantum_numbers_write_array (qn, unit, col_verbose)
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: col_verbose
    integer :: i
    integer :: u
    logical :: col_verb
    u = given_output_unit (unit);  if (u < 0)  return
    col_verb = .false.;  if (present (col_verbose)) col_verb = col_verbose
    write (u, "(A)", advance="no")  "["
    do i = 1, size (qn)
       if (i > 1)  write (u, "(A)", advance="no")  " / "
       if (qn(i)%f%is_defined ()) then
          call qn(i)%f%write (u)
          if (qn(i)%c%is_nonzero () .or. qn(i)%h%is_defined ()) &
               write (u, "(1x)", advance="no")
       end if
       if (col_verb) then
          if (qn(i)%c%is_defined () .or. qn(i)%c%is_ghost ()) then
             call color_write (qn(i)%c, u)
             if (qn(i)%h%is_defined ())  write (u, "(1x)", advance="no")
          end if
       else
          if (qn(i)%c%is_nonzero () .or. qn(i)%c%is_ghost ()) then
             call color_write (qn(i)%c, u)
             if (qn(i)%h%is_defined ())  write (u, "(1x)", advance="no")
          end if
       end if
       if (qn(i)%h%is_defined ()) then
          call qn(i)%h%write (u)
       end if
       if (qn(i)%sub > 0) &
          write (u, "(A,I2)", advance = "no") " SUB = ", qn(i)%sub
    end do
    write (u, "(A)", advance = "no")  "]"
  end subroutine quantum_numbers_write_array

  module subroutine quantum_numbers_write_raw (qn, u)
    class(quantum_numbers_t), intent(in) :: qn
    integer, intent(in) :: u
    call qn%f%write_raw (u)
    call qn%c%write_raw (u)
    call qn%h%write_raw (u)
  end subroutine quantum_numbers_write_raw

  module subroutine quantum_numbers_read_raw (qn, u, iostat)
    class(quantum_numbers_t), intent(out) :: qn
    integer, intent(in) :: u
    integer, intent(out), optional :: iostat
    call qn%f%read_raw (u, iostat=iostat)
    call qn%c%read_raw (u, iostat=iostat)
    call qn%h%read_raw (u, iostat=iostat)
  end subroutine quantum_numbers_read_raw

  impure elemental module function quantum_numbers_get_flavor (qn) result (flv)
    type(flavor_t) :: flv
    class(quantum_numbers_t), intent(in) :: qn
    flv = qn%f
  end function quantum_numbers_get_flavor

  elemental module function quantum_numbers_get_color (qn) result (col)
    type(color_t) :: col
    class(quantum_numbers_t), intent(in) :: qn
    col = qn%c
  end function quantum_numbers_get_color

  elemental module function quantum_numbers_get_helicity (qn) result (hel)
    type(helicity_t) :: hel
    class(quantum_numbers_t), intent(in) :: qn
    hel = qn%h
  end function quantum_numbers_get_helicity

  elemental module function quantum_numbers_get_sub (qn) result (sub)
    integer :: sub
    class(quantum_numbers_t), intent(in) :: qn
    sub = qn%sub
  end function quantum_numbers_get_sub

  elemental module subroutine quantum_numbers_set_color_ghost (qn, ghost)
    class(quantum_numbers_t), intent(inout) :: qn
    logical, intent(in) :: ghost
    call qn%c%set_ghost (ghost)
  end subroutine quantum_numbers_set_color_ghost

  impure elemental module subroutine quantum_numbers_set_model (qn, model)
    class(quantum_numbers_t), intent(inout) :: qn
    class(model_data_t), intent(in), target :: model
    call qn%f%set_model (model)
  end subroutine quantum_numbers_set_model

  elemental module subroutine quantum_numbers_tag_radiated (qn)
    class(quantum_numbers_t), intent(inout) :: qn
    call qn%f%tag_radiated ()
  end subroutine quantum_numbers_tag_radiated

  elemental module subroutine quantum_numbers_tag_hard_process (qn, hard)
    class(quantum_numbers_t), intent(inout) :: qn
    logical, intent(in), optional :: hard
    call qn%f%tag_hard_process (hard)
  end subroutine quantum_numbers_tag_hard_process

  elemental module subroutine quantum_numbers_set_subtraction_index (qn, i)
    class(quantum_numbers_t), intent(inout) :: qn
    integer, intent(in) :: i
    qn%sub = i
  end subroutine quantum_numbers_set_subtraction_index

  elemental module function quantum_numbers_get_subtraction_index &
       (qn) result (sub)
    integer :: sub
    class(quantum_numbers_t), intent(in) :: qn
    sub = qn%sub
  end function quantum_numbers_get_subtraction_index

  elemental module function quantum_numbers_get_color_type (qn) result (color_type)
    integer :: color_type
    class(quantum_numbers_t), intent(in) :: qn
    color_type = qn%f%get_color_type ()
  end function quantum_numbers_get_color_type

  elemental module function quantum_numbers_are_valid (qn) result (valid)
    logical :: valid
    class(quantum_numbers_t), intent(in) :: qn
    valid = qn%f%is_valid ()
  end function quantum_numbers_are_valid

  elemental module function quantum_numbers_are_associated (qn) result (flag)
    logical :: flag
    class(quantum_numbers_t), intent(in) :: qn
    flag = qn%f%is_associated ()
  end function quantum_numbers_are_associated

  elemental module function quantum_numbers_are_diagonal (qn) result (diagonal)
    logical :: diagonal
    class(quantum_numbers_t), intent(in) :: qn
    diagonal = qn%h%is_diagonal () .and. qn%c%is_diagonal ()
  end function quantum_numbers_are_diagonal

  elemental module function quantum_numbers_is_color_ghost (qn) result (ghost)
    logical :: ghost
    class(quantum_numbers_t), intent(in) :: qn
    ghost = qn%c%is_ghost ()
  end function quantum_numbers_is_color_ghost

  elemental module function quantum_numbers_are_hard_process &
       (qn) result (hard_process)
    logical :: hard_process
    class(quantum_numbers_t), intent(in) :: qn
    hard_process = qn%f%is_hard_process ()
  end function quantum_numbers_are_hard_process

  elemental module function quantum_numbers_match (qn1, qn2) result (match)
    logical :: match
    class(quantum_numbers_t), intent(in) :: qn1, qn2
    match = (qn1%f .match. qn2%f) .and. &
         (qn1%c .match. qn2%c) .and. &
         (qn1%h .match. qn2%h)
  end function quantum_numbers_match

  elemental module function quantum_numbers_match_f (qn1, qn2) result (match)
    logical :: match
    class(quantum_numbers_t), intent(in) :: qn1, qn2
    match = (qn1%f .match. qn2%f)
  end function quantum_numbers_match_f

  elemental module function quantum_numbers_match_h (qn1, qn2) result (match)
    logical :: match
    class(quantum_numbers_t), intent(in) :: qn1, qn2
    match = (qn1%h .match. qn2%h)
  end function quantum_numbers_match_h

  elemental module function quantum_numbers_match_fh (qn1, qn2) result (match)
    logical :: match
    class(quantum_numbers_t), intent(in) :: qn1, qn2
    match = (qn1%f .match. qn2%f) .and. &
         (qn1%h .match. qn2%h)
  end function quantum_numbers_match_fh

  elemental module function quantum_numbers_match_hel_diag (qn1, qn2) result (match)
    logical :: match
    class(quantum_numbers_t), intent(in) :: qn1, qn2
    match = (qn1%f .match. qn2%f) .and. &
         (qn1%c .match. qn2%c) .and. &
         (qn1%h .dmatch. qn2%h)
  end function quantum_numbers_match_hel_diag

  elemental module function quantum_numbers_eq_wo_sub (qn1, qn2) result (eq)
    logical :: eq
    type(quantum_numbers_t), intent(in) :: qn1, qn2
    eq = (qn1%f == qn2%f) .and. &
         (qn1%c == qn2%c) .and. &
         (qn1%h == qn2%h)
  end function quantum_numbers_eq_wo_sub

  elemental module function quantum_numbers_eq (qn1, qn2) result (eq)
    logical :: eq
    class(quantum_numbers_t), intent(in) :: qn1, qn2
    eq = (qn1%f == qn2%f) .and. &
         (qn1%c == qn2%c) .and. &
         (qn1%h == qn2%h) .and. &
         (qn1%sub == qn2%sub)
  end function quantum_numbers_eq

  elemental module function quantum_numbers_neq (qn1, qn2) result (neq)
    logical :: neq
    class(quantum_numbers_t), intent(in) :: qn1, qn2
    neq = (qn1%f /= qn2%f) .or. &
         (qn1%c /= qn2%c) .or. &
         (qn1%h /= qn2%h) .or. &
         (qn1%sub /= qn2%sub)
  end function quantum_numbers_neq

  module subroutine quantum_numbers_assign (qn_out, qn_in)
    type(quantum_numbers_t), intent(out) :: qn_out
    type(quantum_numbers_t), intent(in) :: qn_in
    qn_out%f = qn_in%f
    qn_out%c = qn_in%c
    qn_out%h = qn_in%h
    qn_out%sub = qn_in%sub
  end subroutine quantum_numbers_assign

  elemental module function quantum_numbers_are_compatible &
       (qn1, qn2, mask) result (flag)
    logical :: flag
    type(quantum_numbers_t), intent(in) :: qn1, qn2
    type(quantum_numbers_mask_t), intent(in) :: mask
    if (mask%h .or. mask%hd) then
       flag = (qn1%f .match. qn2%f) .and. (qn1%h .match. qn2%h)
    else
       flag = (qn1%f .match. qn2%f)
    end if
    if (mask%c) then
       flag = flag .and. (qn1%c%is_ghost () .eqv. qn2%c%is_ghost ())
    else
       flag = flag .and. &
            .not. (qn1%c%is_ghost () .or. qn2%c%is_ghost ()) .and. &
            (qn1%c == qn2%c)
    end if
  end function quantum_numbers_are_compatible

  elemental module function quantum_numbers_are_physical (qn, mask) result (flag)
    logical :: flag
    type(quantum_numbers_t), intent(in) :: qn
    type(quantum_numbers_mask_t), intent(in) :: mask
    if (mask%c) then
       flag = .true.
    else
       flag = .not. qn%c%is_ghost ()
    end if
  end function quantum_numbers_are_physical

  module subroutine quantum_numbers_canonicalize_color (qn)
    type(quantum_numbers_t), dimension(:), intent(inout) :: qn
    call color_canonicalize (qn%c)
  end subroutine quantum_numbers_canonicalize_color

  module subroutine quantum_numbers_make_color_map (map, qn1, qn2)
    integer, dimension(:,:), intent(out), allocatable :: map
    type(quantum_numbers_t), dimension(:), intent(in) :: qn1, qn2
    call make_color_map (map, qn1%c, qn2%c)
  end subroutine quantum_numbers_make_color_map

  module subroutine quantum_numbers_translate_color0 (qn, map, offset)
    type(quantum_numbers_t), intent(inout) :: qn
    integer, dimension(:,:), intent(in) :: map
    integer, intent(in), optional :: offset
    call color_translate (qn%c, map, offset)
  end subroutine quantum_numbers_translate_color0

  module subroutine quantum_numbers_translate_color1 (qn, map, offset)
    type(quantum_numbers_t), dimension(:), intent(inout) :: qn
    integer, dimension(:,:), intent(in) :: map
    integer, intent(in), optional :: offset
    call color_translate (qn%c, map, offset)
  end subroutine quantum_numbers_translate_color1

  pure module function quantum_numbers_get_max_color_value0 (qn) result (cmax)
    integer :: cmax
    type(quantum_numbers_t), intent(in) :: qn
    cmax = color_get_max_value (qn%c)
  end function quantum_numbers_get_max_color_value0

  pure module function quantum_numbers_get_max_color_value1 (qn) result (cmax)
    integer :: cmax
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    cmax = color_get_max_value (qn%c)
  end function quantum_numbers_get_max_color_value1

  pure module function quantum_numbers_get_max_color_value2 (qn) result (cmax)
    integer :: cmax
    type(quantum_numbers_t), dimension(:,:), intent(in) :: qn
    cmax = color_get_max_value (qn%c)
  end function quantum_numbers_get_max_color_value2

  elemental module subroutine quantum_numbers_add_color_offset (qn, offset)
    class(quantum_numbers_t), intent(inout) :: qn
    integer, intent(in) :: offset
    call qn%c%add_offset (offset)
  end subroutine quantum_numbers_add_color_offset

  module subroutine quantum_number_array_make_color_contractions (qn_in, qn_out)
    type(quantum_numbers_t), dimension(:), intent(in) :: qn_in
    type(quantum_numbers_t), dimension(:,:), intent(out), allocatable :: qn_out
    type(color_t), dimension(:,:), allocatable :: col
    integer :: i
    call color_array_make_contractions (qn_in%c, col)
    allocate (qn_out (size (col, 1), size (col, 2)))
    do i = 1, size (qn_out, 2)
       qn_out(:,i)%f = qn_in%f
       qn_out(:,i)%c = col(:,i)
       qn_out(:,i)%h = qn_in%h
    end do
  end subroutine quantum_number_array_make_color_contractions

  elemental module subroutine quantum_numbers_invert_color (qn)
    class(quantum_numbers_t), intent(inout) :: qn
    call qn%c%invert ()
  end subroutine quantum_numbers_invert_color

  elemental module subroutine quantum_numbers_flip_helicity (qn)
    class(quantum_numbers_t), intent(inout) :: qn
    call qn%h%flip ()
  end subroutine quantum_numbers_flip_helicity

  module function merge_quantum_numbers0 (qn1, qn2) result (qn3)
    type(quantum_numbers_t) :: qn3
    type(quantum_numbers_t), intent(in) :: qn1, qn2
    qn3%f = qn1%f .merge. qn2%f
    qn3%c = qn1%c .merge. qn2%c
    qn3%h = qn1%h .merge. qn2%h
    qn3%sub = merge_subtraction_index (qn1%sub, qn2%sub)
  end function merge_quantum_numbers0

  module function merge_quantum_numbers1 (qn1, qn2) result (qn3)
    type(quantum_numbers_t), dimension(:), intent(in) :: qn1, qn2
    type(quantum_numbers_t), dimension(size(qn1)) :: qn3
    qn3%f = qn1%f .merge. qn2%f
    qn3%c = qn1%c .merge. qn2%c
    qn3%h = qn1%h .merge. qn2%h
    qn3%sub = merge_subtraction_index (qn1%sub, qn2%sub)
  end function merge_quantum_numbers1

  elemental function merge_subtraction_index (sub1, sub2) result (sub3)
    integer :: sub3
    integer, intent(in) :: sub1, sub2
    if (sub1 > 0 .and. sub2 > 0) then
       if (sub1 == sub2) then
          sub3 = sub1
       else
          sub3 = 0
       end if
    else if (sub1 > 0) then
       sub3 = sub1
    else if (sub2 > 0) then
       sub3 = sub2
    else
       sub3 = 0
    end if
  end function merge_subtraction_index

  elemental module function quantum_numbers_mask &
       (mask_f, mask_c, mask_h, mask_cg, mask_hd) result (mask)
    type(quantum_numbers_mask_t) :: mask
    logical, intent(in) :: mask_f, mask_c, mask_h
    logical, intent(in), optional :: mask_cg
    logical, intent(in), optional :: mask_hd
    call quantum_numbers_mask_init &
         (mask, mask_f, mask_c, mask_h, mask_cg, mask_hd)
  end function quantum_numbers_mask

  elemental module subroutine quantum_numbers_mask_init &
       (mask, mask_f, mask_c, mask_h, mask_cg, mask_hd)
    class(quantum_numbers_mask_t), intent(inout) :: mask
    logical, intent(in) :: mask_f, mask_c, mask_h
    logical, intent(in), optional :: mask_cg, mask_hd
    mask%f = mask_f
    mask%c = mask_c
    mask%h = mask_h
    mask%cg = .false.
    if (present (mask_cg)) then
       if (mask%c)  mask%cg = mask_cg
    else
       mask%cg = mask_c
    end if
    mask%hd = .false.
    if (present (mask_hd)) then
       if (.not. mask%h)  mask%hd = mask_hd
    end if
  end subroutine quantum_numbers_mask_init

  module subroutine quantum_numbers_mask_write_single (mask, unit)
    class(quantum_numbers_mask_t), intent(in) :: mask
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(A)", advance="no") "["
    write (u, "(L1)", advance="no")  mask%f
    write (u, "(L1)", advance="no")  mask%c
    if (.not.mask%cg)  write (u, "('g')", advance="no")
    write (u, "(L1)", advance="no")  mask%h
    if (mask%hd)  write (u, "('d')", advance="no")
    write (u, "(A)", advance="no") "]"
  end subroutine quantum_numbers_mask_write_single

  module subroutine quantum_numbers_mask_write_array (mask, unit)
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: mask
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(A)", advance="no") "["
    do i = 1, size (mask)
       if (i > 1)  write (u, "(A)", advance="no")  "/"
       write (u, "(L1)", advance="no")  mask(i)%f
       write (u, "(L1)", advance="no")  mask(i)%c
       if (.not.mask(i)%cg)  write (u, "('g')", advance="no")
       write (u, "(L1)", advance="no")  mask(i)%h
       if (mask(i)%hd)  write (u, "('d')", advance="no")
    end do
    write (u, "(A)", advance="no") "]"
  end subroutine quantum_numbers_mask_write_array

  elemental module subroutine quantum_numbers_mask_set_flavor (mask, mask_f)
    class(quantum_numbers_mask_t), intent(inout) :: mask
    logical, intent(in) :: mask_f
    mask%f = mask_f
  end subroutine quantum_numbers_mask_set_flavor

  elemental module subroutine quantum_numbers_mask_set_color (mask, mask_c, mask_cg)
    class(quantum_numbers_mask_t), intent(inout) :: mask
    logical, intent(in) :: mask_c
    logical, intent(in), optional :: mask_cg
    mask%c = mask_c
    if (present (mask_cg)) then
       if (mask%c)  mask%cg = mask_cg
    else
       mask%cg = mask_c
    end if
  end subroutine quantum_numbers_mask_set_color

  elemental module subroutine quantum_numbers_mask_set_helicity (mask, mask_h, mask_hd)
    class(quantum_numbers_mask_t), intent(inout) :: mask
    logical, intent(in) :: mask_h
    logical, intent(in), optional :: mask_hd
    mask%h = mask_h
    if (present (mask_hd)) then
       if (.not. mask%h)  mask%hd = mask_hd
    end if
  end subroutine quantum_numbers_mask_set_helicity

  elemental module subroutine quantum_numbers_mask_set_sub (mask, sub)
    class(quantum_numbers_mask_t), intent(inout) :: mask
    integer, intent(in) :: sub
    mask%sub = sub
  end subroutine quantum_numbers_mask_set_sub

  elemental module subroutine quantum_numbers_mask_assign &
       (mask, mask_in, flavor, color, helicity)
    class(quantum_numbers_mask_t), intent(inout) :: mask
    class(quantum_numbers_mask_t), intent(in) :: mask_in
    logical, intent(in), optional :: flavor, color, helicity
    if (present (flavor)) then
       if (flavor) then
          mask%f = mask_in%f
       end if
    end if
    if (present (color)) then
       if (color) then
          mask%c = mask_in%c
          mask%cg = mask_in%cg
       end if
    end if
    if (present (helicity)) then
       if (helicity) then
          mask%h = mask_in%h
          mask%hd = mask_in%hd
       end if
    end if
  end subroutine quantum_numbers_mask_assign

  module function quantum_numbers_mask_any (mask) result (match)
    logical :: match
    type(quantum_numbers_mask_t), intent(in) :: mask
    match = mask%f .or. mask%c .or. mask%h .or. mask%hd
  end function quantum_numbers_mask_any

  elemental module function quantum_numbers_mask_or (mask1, mask2) result (mask)
    type(quantum_numbers_mask_t) :: mask
    class(quantum_numbers_mask_t), intent(in) :: mask1, mask2
    mask%f = mask1%f .or. mask2%f
    mask%c = mask1%c .or. mask2%c
    if (mask%c)  mask%cg = mask1%cg .or. mask2%cg
    mask%h = mask1%h .or. mask2%h
    if (.not. mask%h)  mask%hd = mask1%hd .or. mask2%hd
  end function quantum_numbers_mask_or

  elemental module function quantum_numbers_mask_eqv (mask1, mask2) result (eqv)
    logical :: eqv
    class(quantum_numbers_mask_t), intent(in) :: mask1, mask2
    eqv = (mask1%f .eqv. mask2%f) .and. &
         (mask1%c .eqv. mask2%c) .and. &
         (mask1%cg .eqv. mask2%cg) .and. &
         (mask1%h .eqv. mask2%h) .and. &
         (mask1%hd .eqv. mask2%hd)
  end function quantum_numbers_mask_eqv

  elemental module function quantum_numbers_mask_neqv (mask1, mask2) result (neqv)
    logical :: neqv
    class(quantum_numbers_mask_t), intent(in) :: mask1, mask2
    neqv = (mask1%f .neqv. mask2%f) .or. &
         (mask1%c .neqv. mask2%c) .or. &
         (mask1%cg .neqv. mask2%cg) .or. &
         (mask1%h .neqv. mask2%h) .or. &
         (mask1%hd .neqv. mask2%hd)
  end function quantum_numbers_mask_neqv

  elemental module subroutine quantum_numbers_undefine (qn, mask)
    class(quantum_numbers_t), intent(inout) :: qn
    type(quantum_numbers_mask_t), intent(in) :: mask
    if (mask%f)  call qn%f%undefine ()
    if (mask%c)  call qn%c%undefine (undefine_ghost = mask%cg)
    if (mask%h) then
       call qn%h%undefine ()
    else if (mask%hd) then
       if (.not. qn%h%is_diagonal ()) then
          call qn%h%diagonalize ()
       end if
    end if
    if (mask%sub > 0)  qn%sub = 0
  end subroutine quantum_numbers_undefine

  module function quantum_numbers_undefined0 (qn, mask) result (qn_new)
    class(quantum_numbers_t), intent(in) :: qn
    type(quantum_numbers_mask_t), intent(in) :: mask
    type(quantum_numbers_t) :: qn_new
    select type (qn)
    type is (quantum_numbers_t);  qn_new = qn
    end select
    call quantum_numbers_undefine (qn_new, mask)
  end function quantum_numbers_undefined0

  module function quantum_numbers_undefined1 (qn, mask) result (qn_new)
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    type(quantum_numbers_mask_t), intent(in) :: mask
    type(quantum_numbers_t), dimension(size(qn)) :: qn_new
    qn_new = qn
    call quantum_numbers_undefine (qn_new, mask)
  end function quantum_numbers_undefined1

  module function quantum_numbers_undefined11 (qn, mask) result (qn_new)
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: mask
    type(quantum_numbers_t), dimension(size(qn)) :: qn_new
    qn_new = qn
    call quantum_numbers_undefine (qn_new, mask)
  end function quantum_numbers_undefined11

  elemental module function quantum_numbers_are_redundant (qn, mask) &
       result (redundant)
    logical :: redundant
    class(quantum_numbers_t), intent(in) :: qn
    type(quantum_numbers_mask_t), intent(in) :: mask
    redundant = .false.
    if (mask%f) then
       redundant = qn%f%is_defined ()
    end if
    if (mask%c) then
       redundant = qn%c%is_defined ()
    end if
    if (mask%h) then
       redundant = qn%h%is_defined ()
    else if (mask%hd) then
       redundant = .not. qn%h%is_diagonal ()
    end if
    if (mask%sub > 0) redundant = qn%sub >= mask%sub
  end function quantum_numbers_are_redundant

  elemental module function quantum_numbers_mask_diagonal_helicity (mask) &
       result (flag)
    logical :: flag
    class(quantum_numbers_mask_t), intent(in) :: mask
    flag = mask%h .or. mask%hd
  end function quantum_numbers_mask_diagonal_helicity


end submodule quantum_numbers_s

