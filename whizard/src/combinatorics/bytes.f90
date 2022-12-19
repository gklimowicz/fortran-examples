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

module bytes

  use kinds, only: i8, i32, i64

  implicit none
  private

  public :: byte_t
  public :: byte_zero
  public :: assignment(=)
  public :: byte_write
  public :: word32_t
  public :: word32_empty, word32_filled, word32_fill
  public :: word32_append_byte
  public :: byte_from_word32
  public :: word32_write
  public :: not, ior, ieor, iand, ishft, ishftc
  public :: operator(+)
  public :: word64_t
  public :: byte_from_word64, word32_from_word64
  public :: word64_write

  type :: byte_t
     private
     integer(i8) :: i
  end type byte_t

  type :: word32_t
     private
     integer(i32) :: i
     integer :: fill = 0
  end type word32_t

  type :: word64_t
     private
     integer(i64) :: i
  end type word64_t


  type(byte_t), parameter :: byte_zero = byte_t (0_i8)


  interface assignment(=)
     module procedure set_byte_from_i8
  end interface
  interface byte_write
     module procedure byte_write_unit, byte_write_string
  end interface
  interface assignment(=)
     module procedure word32_set_from_i32
     module procedure word32_set_from_byte
  end interface
  interface assignment(=)
     module procedure i32_from_word32
  end interface
  interface word32_write
     module procedure word32_write_unit
  end interface
  interface not
     module procedure word_not
  end interface
  interface ior
     module procedure word_or
  end interface
  interface ieor
     module procedure word_eor
  end interface
  interface iand
     module procedure word_and
  end interface
  interface ishft
     module procedure word_shft
  end interface
  interface ishftc
     module procedure word_shftc
  end interface
  interface operator(+)
     module procedure word_add
     module procedure word_add_i8
     module procedure word_add_i32
  end interface
  interface assignment(=)
     module procedure word64_set_from_i64
     module procedure word64_set_from_word32
  end interface
  interface word64_write
     module procedure word64_write_unit
  end interface

  interface
    module subroutine set_byte_from_i8 (b, i)
      type(byte_t), intent(out) :: b
      integer(i8), intent(in) :: i
    end subroutine set_byte_from_i8
    module subroutine byte_write_unit (b, unit, decimal, newline)
      type(byte_t), intent(in), optional :: b
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: decimal, newline
    end subroutine byte_write_unit
    module subroutine byte_write_string (b, s)
      type(byte_t), intent(in) :: b
      character(len=2), intent(inout) :: s
    end subroutine byte_write_string
    module subroutine word32_set_from_i32 (w, i)
      type(word32_t), intent(out) :: w
      integer(i32), intent(in) :: i
    end subroutine word32_set_from_i32
    module subroutine i32_from_word32 (i, w)
      integer(i32), intent(out) :: i
      type(word32_t), intent(in) :: w
    end subroutine i32_from_word32
    module subroutine word32_set_from_byte (w, b)
      type(word32_t), intent(out) :: w
      type(byte_t), intent(in) :: b
    end subroutine word32_set_from_byte
    module function word32_empty (w)
      type(word32_t), intent(in) :: w
      logical :: word32_empty
    end function word32_empty
    module function word32_filled (w)
      type(word32_t), intent(in) :: w
      logical :: word32_filled
    end function word32_filled
    module function word32_fill (w)
      type(word32_t), intent(in) :: w
      integer :: word32_fill
    end function word32_fill
    module subroutine word32_append_byte (w, b)
      type(word32_t), intent(inout) :: w
      type(byte_t), intent(in) :: b
    end subroutine word32_append_byte
    module function byte_from_word32 (w, i) result (b)
      type(word32_t), intent(in) :: w
      integer, intent(in) :: i
      type(byte_t) :: b
    end function byte_from_word32
    module subroutine word32_write_unit (w, unit, bytes, decimal, newline)
      type(word32_t), intent(in) :: w
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: bytes, decimal, newline
    end subroutine word32_write_unit
    module function word_not (w1) result (w2)
      type(word32_t), intent(in) :: w1
      type(word32_t) :: w2
    end function word_not
    module function word_or (w1, w2) result (w3)
      type(word32_t), intent(in) :: w1, w2
      type(word32_t) :: w3
    end function word_or
    module function word_eor (w1, w2) result (w3)
      type(word32_t), intent(in) :: w1, w2
      type(word32_t) :: w3
    end function word_eor
    module function word_and (w1, w2) result (w3)
      type(word32_t), intent(in) :: w1, w2
      type(word32_t) :: w3
    end function word_and
    module function word_shft (w1, s) result (w2)
      type(word32_t), intent(in) :: w1
      integer, intent(in) :: s
      type(word32_t) :: w2
    end function word_shft
    module function word_shftc (w1, s) result (w2)
      type(word32_t), intent(in) :: w1
      integer, intent(in) :: s
      type(word32_t) :: w2
    end function word_shftc
  module function word_add (w1, w2) result (w3)
    type(word32_t), intent(in) :: w1, w2
    type(word32_t) :: w3
  end function word_add
  module function word_add_i8 (w1, i) result (w3)
    type(word32_t), intent(in) :: w1
    integer(i8), intent(in) :: i
    type(word32_t) :: w3
  end function word_add_i8
  module function word_add_i32 (w1, i) result (w3)
    type(word32_t), intent(in) :: w1
    integer(i32), intent(in) :: i
    type(word32_t) :: w3
  end function word_add_i32
    module subroutine word64_set_from_i64 (ww, i)
      type(word64_t), intent(out) :: ww
      integer(i64), intent(in) :: i
    end subroutine word64_set_from_i64
    module subroutine word64_set_from_word32 (ww, w)
      type(word64_t), intent(out) :: ww
      type(word32_t), intent(in) :: w
    end subroutine word64_set_from_word32
    module function byte_from_word64 (ww, i) result (b)
      type(word64_t), intent(in) :: ww
      integer, intent(in) :: i
      type(byte_t) :: b
    end function byte_from_word64
    module function word32_from_word64 (ww, i) result (w)
      type(word64_t), intent(in) :: ww
      integer, intent(in) :: i
      type(word32_t) :: w
    end function word32_from_word64
    module subroutine word64_write_unit (ww, unit, words, bytes, decimal, newline)
      type(word64_t), intent(in) :: ww
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: words, bytes, decimal, newline
    end subroutine word64_write_unit
  end interface

end module bytes
