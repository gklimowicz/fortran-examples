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

submodule (bytes) bytes_s

  use io_units

  implicit none

contains

  module subroutine set_byte_from_i8 (b, i)
    type(byte_t), intent(out) :: b
    integer(i8), intent(in) :: i
    b%i = i
  end subroutine set_byte_from_i8

  module subroutine byte_write_unit (b, unit, decimal, newline)
    type(byte_t), intent(in), optional :: b
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: decimal, newline
    logical :: dc, nl
    type(word32_t) :: w
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    dc = .false.;  if (present (decimal))  dc = decimal
    nl = .false.;  if (present (newline))  nl = newline
    if (dc) then
       w = b
       write (u, '(I3)', advance='no')  w%i
    else
       write (u, '(z2.2)', advance='no')  b%i
    end if
    if (nl) write (u, *)
  end subroutine byte_write_unit

  module subroutine byte_write_string (b, s)
    type(byte_t), intent(in) :: b
    character(len=2), intent(inout) :: s
    write (s, '(z2.2)')  b%i
  end subroutine byte_write_string

  module subroutine word32_set_from_i32 (w, i)
    type(word32_t), intent(out) :: w
    integer(i32), intent(in) :: i
    w%i = i
    w%fill = 32
  end subroutine word32_set_from_i32

  module subroutine i32_from_word32 (i, w)
    integer(i32), intent(out) :: i
    type(word32_t), intent(in) :: w
    i = w%i
  end subroutine i32_from_word32

  module subroutine word32_set_from_byte (w, b)
    type(word32_t), intent(out) :: w
    type(byte_t), intent(in) :: b
    if (b%i >= 0_i8) then
       w%i = b%i
    else
       w%i = 2_i32*(huge(0_i8)+1_i32) + b%i
    end if
    w%fill = 32
  end subroutine word32_set_from_byte

  module function word32_empty (w)
    type(word32_t), intent(in) :: w
    logical :: word32_empty
    word32_empty = (w%fill == 0)
  end function word32_empty

  module function word32_filled (w)
    type(word32_t), intent(in) :: w
    logical :: word32_filled
    word32_filled = (w%fill == 32)
  end function word32_filled

  module function word32_fill (w)
    type(word32_t), intent(in) :: w
    integer :: word32_fill
    word32_fill = w%fill
  end function word32_fill

  module subroutine word32_append_byte (w, b)
    type(word32_t), intent(inout) :: w
    type(byte_t), intent(in) :: b
    type(word32_t) :: w1
    if (.not. word32_filled (w)) then
       w1 = b
       call mvbits (w1%i, 0, 8, w%i, w%fill)
       w%fill = w%fill + 8
    end if
  end subroutine word32_append_byte

  module function byte_from_word32 (w, i) result (b)
    type(word32_t), intent(in) :: w
    integer, intent(in) :: i
    type(byte_t) :: b
    integer(i32) :: j
    j = 0
    if (i >= 0 .and. i*8 < w%fill) then
       call mvbits (w%i, i*8, 8, j, 0)
    end if
    b%i = int (ibclr (j, 7), kind=i8)
    if (btest (j, 7))  b%i = ibset (b%i, 7)
  end function byte_from_word32

  module subroutine word32_write_unit (w, unit, bytes, decimal, newline)
    type(word32_t), intent(in) :: w
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: bytes, decimal, newline
    logical :: dc, by, nl
    type(word64_t) :: ww
    integer :: i, u
    u = given_output_unit (unit);  if (u < 0)  return
    by = .false.;  if (present (bytes))   by = bytes
    dc = .false.;  if (present (decimal)) dc = decimal
    nl = .false.;  if (present (newline)) nl = newline
    if (by) then
       do i = 0, 3
          if (i>0)  write (u, '(1x)', advance='no')
          if (8*i < w%fill) then
             call byte_write (byte_from_word32 (w, i), unit, decimal=decimal)
          else if (dc) then
             write (u, '(3x)', advance='no')
          else
             write (u, '(2x)', advance='no')
          end if
       end do
    else if (dc) then
       ww = w
       write (u, '(I10)', advance='no') ww%i
    else
       select case (w%fill)
       case ( 0)
       case ( 8);  write (6, '(1x,z8.2)', advance='no') ibits (w%i, 0, 8)
       case (16);  write (6, '(1x,z8.4)', advance='no') ibits (w%i, 0,16)
       case (24);  write (6, '(1x,z8.6)', advance='no') ibits (w%i, 0,24)
       case (32);  write (6, '(1x,z8.8)', advance='no') ibits (w%i, 0,32)
       end select
    end if
    if (nl) write (u, *)
  end subroutine word32_write_unit

  module function word_not (w1) result (w2)
    type(word32_t), intent(in) :: w1
    type(word32_t) :: w2
    w2 = not (w1%i)
  end function word_not

  module function word_or (w1, w2) result (w3)
    type(word32_t), intent(in) :: w1, w2
    type(word32_t) :: w3
    w3 = ior (w1%i, w2%i)
  end function word_or

  module function word_eor (w1, w2) result (w3)
    type(word32_t), intent(in) :: w1, w2
    type(word32_t) :: w3
    w3 = ieor (w1%i, w2%i)
  end function word_eor

  module function word_and (w1, w2) result (w3)
    type(word32_t), intent(in) :: w1, w2
    type(word32_t) :: w3
    w3 = iand (w1%i, w2%i)
  end function word_and

  module function word_shft (w1, s) result (w2)
    type(word32_t), intent(in) :: w1
    integer, intent(in) :: s
    type(word32_t) :: w2
    w2 = ishft (w1%i, s)
  end function word_shft

  module function word_shftc (w1, s) result (w2)
    type(word32_t), intent(in) :: w1
    integer, intent(in) :: s
    type(word32_t) :: w2
    w2 = ishftc (w1%i, s, 32)
  end function word_shftc

  module function word_add (w1, w2) result (w3)
    type(word32_t), intent(in) :: w1, w2
    type(word32_t) :: w3
    integer(i64) :: j
    j = int (ibclr (w1%i, 31), i64) + int (ibclr (w2%i, 31), i64)
    w3 = int (ibclr (j, 31), kind=i32)
    if (btest (j, 31)) then
       if (btest (w1%i, 31) .eqv. btest (w2%i, 31))  w3 = ibset (w3%i, 31)
    else
       if (btest (w1%i, 31) .neqv. btest (w2%i, 31))  w3 = ibset (w3%i, 31)
    end if
  end function word_add

  module function word_add_i8 (w1, i) result (w3)
    type(word32_t), intent(in) :: w1
    integer(i8), intent(in) :: i
    type(word32_t) :: w3
    integer(i64) :: j
    j = int (ibclr (w1%i, 31), i64) + int (ibclr (i, 7), i64)
    if (btest (i, 7))  j = j + 128
    w3 = int (ibclr (j, 31), kind=i32)
    if (btest (j, 31) .neqv. btest (w1%i, 31))  w3 = ibset (w3%i, 31)
  end function word_add_i8

  module function word_add_i32 (w1, i) result (w3)
    type(word32_t), intent(in) :: w1
    integer(i32), intent(in) :: i
    type(word32_t) :: w3
    integer(i64) :: j
    j = int (ibclr (w1%i, 31), i64) + int (ibclr (i, 31), i64)
    w3 = int (ibclr (j, 31), kind=i32)
    if (btest (j, 31)) then
       if (btest (w1%i, 31) .eqv. btest (i, 31))  w3 = ibset (w3%i, 31)
    else
       if (btest (w1%i, 31) .neqv. btest (i, 31))  w3 = ibset (w3%i, 31)
    end if
  end function word_add_i32

  module subroutine word64_set_from_i64 (ww, i)
    type(word64_t), intent(out) :: ww
    integer(i64), intent(in) :: i
    ww%i = i
  end subroutine word64_set_from_i64

  module subroutine word64_set_from_word32 (ww, w)
    type(word64_t), intent(out) :: ww
    type(word32_t), intent(in) :: w
    if (w%i >= 0_i32) then
       ww%i = w%i
    else
       ww%i = 2_i64*(huge(0_i32)+1_i64) + w%i
    end if
  end subroutine word64_set_from_word32

  module function byte_from_word64 (ww, i) result (b)
    type(word64_t), intent(in) :: ww
    integer, intent(in) :: i
    type(byte_t) :: b
    integer(i64) :: j
    j = 0
    if (i >= 0 .and. i*8 < 64) then
       call mvbits (ww%i, i*8, 8, j, 0)
    end if
    b%i = int (ibclr (j, 7), kind=i8)
    if (btest (j, 7))  b%i = ibset (b%i, 7)
  end function byte_from_word64

  module function word32_from_word64 (ww, i) result (w)
    type(word64_t), intent(in) :: ww
    integer, intent(in) :: i
    type(word32_t) :: w
    integer(i64) :: j
    j = 0
    select case (i)
    case (0);  call mvbits (ww%i,  0, 32, j, 0)
    case (1);  call mvbits (ww%i, 32, 32, j, 0)
    end select
    w = int (ibclr (j, 31), kind=i32)
    if (btest (j, 31))  w = ibset (w%i, 31)
  end function word32_from_word64

  module subroutine word64_write_unit (ww, unit, words, bytes, decimal, newline)
    type(word64_t), intent(in) :: ww
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: words, bytes, decimal, newline
    logical :: wo, by, dc, nl
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    wo = .false.;  if (present (words))    wo = words
    by = .false.;  if (present (bytes))    by = bytes
    dc = .false.;  if (present (decimal))  dc = decimal
    nl = .false.;  if (present (newline))  nl = newline
    if (wo .or. by) then
       call word32_write_unit (word32_from_word64 (ww, 0), unit, by, dc)
       write (u, '(2x)', advance='no')
       call word32_write_unit (word32_from_word64 (ww, 1), unit, by, dc)
    else if (dc) then
       write (u, '(I19)', advance='no') ww%i
    else
       write (u, '(Z16)', advance='no') ww%i
    end if
    if (nl) write (u, *)
  end subroutine word64_write_unit


end submodule bytes_s

