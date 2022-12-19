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

submodule (md5) md5_s

  use kinds, only: i8, i32, i64
  use system_defs, only: BUFFER_SIZE
  use system_defs, only: LF, EOR, EOF
  use io_units
  use diagnostics

  implicit none

  type(word32_t), dimension(64), save :: t
  logical, save :: table_initialized = .false.

contains

  function block_is_empty (b)
    type(block_t), intent(in) :: b
    logical :: block_is_empty
    block_is_empty = (b%fill == 0 .and. word32_empty (b%w(0)))
  end function block_is_empty

  function block_is_filled (b)
    type(block_t), intent(in) :: b
    logical :: block_is_filled
    block_is_filled = (b%fill == 64)
  end function block_is_filled

  subroutine block_append_byte (bl, by)
    type(block_t), intent(inout) :: bl
    type(byte_t), intent(in) :: by
    if (.not. block_is_filled (bl)) then
       call word32_append_byte (bl%w(bl%fill/4), by)
       bl%fill = bl%fill + 1
    end if
  end subroutine block_append_byte

  module subroutine block_write_unit (b, unit, bytes, decimal)
    type(block_t), intent(in) :: b
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: bytes, decimal
    logical :: by, dc
    integer :: i, u
    u = given_output_unit (unit);  if (u < 0)  return
    by = .false.;  if (present (bytes))    by = bytes
    dc = .false.;  if (present (decimal))  dc = decimal
    do i = 0, b%fill/4 - 1
       call newline_or_blank (u, i, by, dc)
       call word32_write (b%w(i), unit, bytes, decimal)
    end do
    if (.not. block_is_filled (b)) then
       i = b%fill/4
       if (.not. word32_empty (b%w(i))) then
          call newline_or_blank (u, i, by, dc)
          call word32_write (b%w(i), unit, bytes, decimal)
       end if
    end if
    write (u, *)
  contains
    subroutine newline_or_blank (u, i, bytes, decimal)
      integer, intent(in) :: u, i
      logical, intent(in) :: bytes, decimal
      if (decimal) then
         select case (i)
         case (0)
         case (2,4,6,8,10,12,14);  write (u, *)
         case default
            write (u, '(2x)', advance='no')
         end select
      else if (bytes) then
         select case (i)
         case (0)
         case (4,8,12);  write (u, *)
         case default
            write (u, '(2x)', advance='no')
         end select
      else
         if (i == 8)  write (u, *)
      end if
    end subroutine newline_or_blank
  end subroutine block_write_unit

  subroutine message_clear (m)
    type(message_t), intent(inout) :: m
    type(block_t), pointer :: b
    nullify (m%last)
    do
       b => m%first
       if (.not.(associated (b))) exit
       m%first => b%next
       deallocate (b)
    end do
    m%n_blocks = 0
  end subroutine message_clear

  subroutine message_append_new_block (m)
    type(message_t), intent(inout) :: m
    if (associated (m%last)) then
       allocate (m%last%next)
       m%last => m%last%next
       m%n_blocks = m%n_blocks + 1
    else
       allocate (m%first)
       m%last => m%first
       m%n_blocks = 1
    end if
  end subroutine message_append_new_block

  subroutine message_init (m)
    type(message_t), intent(inout) :: m
    call message_clear (m)
    call message_append_new_block (m)
  end subroutine message_init

  subroutine message_append_byte (m, b)
    type(message_t), intent(inout) :: m
    type(byte_t), intent(in) :: b
    if (.not. associated (m%last)) then
       call message_init (m)
    else if (block_is_filled (m%last)) then
       call message_append_new_block (m)
    end if
    call block_append_byte (m%last, b)
  end subroutine message_append_byte

  subroutine message_pad_zero (m, i)
    type(message_t), intent(inout) :: m
    integer, intent(in) :: i
    type(block_t), pointer :: b
    integer :: j
    if (associated (m%last)) then
       b => m%last
       if (b%fill > i) then
          do j = b%fill + 1, 64 + i
             call message_append_byte (m, byte_zero)
          end do
       else
          do j = b%fill + 1, i
             call message_append_byte (m, byte_zero)
          end do
       end if
    end if
  end subroutine message_pad_zero

  function message_bits (m) result (length)
    type(message_t), intent(in) :: m
    type(word64_t) :: length
    type(block_t), pointer :: b
    integer(i64) :: n_blocks_filled, n_bytes_extra
    if (m%n_blocks > 0) then
       b => m%last
       if (block_is_filled (b)) then
          n_blocks_filled = m%n_blocks
          n_bytes_extra = 0
       else
          n_blocks_filled = m%n_blocks - 1
          n_bytes_extra = b%fill
       end if
       length = n_blocks_filled * 512 + n_bytes_extra * 8
    else
       length = 0_i64
    end if
  end function message_bits

  subroutine message_append_string (m, s)
    type(message_t), intent(inout) :: m
    character(len=*), intent(in) :: s
    integer(i64) :: i, n_bytes
    integer(i8), dimension(:), allocatable :: buffer
    integer(i8), dimension(1) :: mold
    type(byte_t) :: b
    n_bytes = size (transfer (s, mold))
    allocate (buffer (n_bytes))
    buffer = transfer (s, mold)
    do i = 1, size (buffer)
       b = buffer(i)
       call message_append_byte (m, b)
    end do
    deallocate (buffer)
  end subroutine message_append_string

  subroutine message_append_i32 (m, x)
    type(message_t), intent(inout) :: m
    integer(i32), intent(in) :: x
    integer(i8), dimension(4) :: buffer
    type(byte_t) :: b
    integer :: i
    buffer = transfer (x, buffer, size(buffer))
    do i = 1, size (buffer)
       b = buffer(i)
       call message_append_byte (m, b)
    end do
  end subroutine message_append_i32

  subroutine message_append_from_unit (m, u, iostat)
    type(message_t), intent(inout) :: m
    integer, intent(in) :: u
    integer, intent(out) :: iostat
    character(len=BUFFER_SIZE) :: buffer
    read (u, *, iostat=iostat) buffer
    call message_append_string (m, trim (buffer))
    call message_append_string (m, LF)
  end subroutine message_append_from_unit

  subroutine message_read_from_file (m, f)
    type(message_t), intent(inout) :: m
    character(len=*), intent(in) :: f
    integer :: u, iostat
    u = free_unit ()
    open (file=f, unit=u, action='read')
    do
       call message_append_from_unit (m, u, iostat=iostat)
       if (iostat < 0) exit
    end do
    close (u)
  end subroutine message_read_from_file

  module subroutine message_write_unit (m, unit, bytes, decimal)
    type(message_t), intent(in) :: m
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: bytes, decimal
    type(block_t), pointer :: b
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    b => m%first
    if (associated (b)) then
       do
          call block_write_unit (b, unit, bytes, decimal)
          b => b%next
          if (.not. associated (b))  exit
          write (u, *)
       end do
    end if
  end subroutine message_write_unit

  function ff (x, y, z)
    type(word32_t), intent(in) :: x, y, z
    type(word32_t) :: ff
    ff = ior (iand (x, y), iand (not (x), z))
  end function ff

  function fg (x, y, z)
    type(word32_t), intent(in) :: x, y, z
    type(word32_t) :: fg
    fg = ior (iand (x, z), iand (y, not (z)))
  end function fg

  function fh (x, y, z)
    type(word32_t), intent(in) :: x, y, z
    type(word32_t) :: fh
    fh = ieor (ieor (x, y), z)
  end function fh

  function fi (x, y, z)
    type(word32_t), intent(in) :: x, y, z
    type(word32_t) :: fi
    fi = ieor (y, ior (x, not (z)))
  end function fi

  subroutine table_init
    type(word64_t) :: ww
    integer :: i
    if (.not.table_initialized) then
       do i = 1, 64
          ww = int (4294967296d0 * abs (sin (i * 1d0)), kind=i64)
          t(i) = word32_from_word64 (ww, 0)
       end do
       table_initialized = .true.
    end if
  end subroutine table_init

  function digest_string (aa) result (s)
    type(word32_t), dimension (0:3), intent(in) :: aa
    character(len=32) :: s
    integer :: i, j
    do i = 0, 3
       do j = 0, 3
          call byte_write (byte_from_word32 (aa(i), j), s(i*8+j*2+1:i*8+j*2+2))
       end do
    end do
  end function digest_string

  subroutine message_pad (m)
    type(message_t), intent(inout) :: m
    type(word64_t) :: length
    integer(i8), parameter :: ipad = -128  ! z'80'
    type(byte_t) :: b
    integer :: i
    length = message_bits (m)
    b = ipad
    call message_append_byte (m, b)
    call message_pad_zero (m, 56)
    do i = 0, 7
       call message_append_byte (m, byte_from_word64 (length, i))
    end do
  end subroutine message_pad

  subroutine message_digest (m, s)
    type(message_t), intent(in) :: m
    character(len=32), intent(out) :: s
    integer(i32), parameter :: ia =  1732584193  ! z'67452301'
    integer(i32), parameter :: ib =  -271733879  ! z'efcdab89'
    integer(i32), parameter :: ic = -1732584194  ! z'98badcfe'
    integer(i32), parameter :: id =   271733878  ! z'10325476'
    type(word32_t) :: a, b, c, d
    type(word32_t) :: aa, bb, cc, dd
    type(word32_t), dimension(0:15) :: x
    type(block_t), pointer :: bl
    call table_init
    a = ia;  b = ib;  c = ic;  d = id
    bl => m%first
    do
       if (.not.associated (bl)) exit
       x = bl%w
       aa = a;  bb = b;  cc = c;  dd = d
       call transform (ff, a, b, c, d,  0,  7,  1)
       call transform (ff, d, a, b, c,  1, 12,  2)
       call transform (ff, c, d, a, b,  2, 17,  3)
       call transform (ff, b, c, d, a,  3, 22,  4)
       call transform (ff, a, b, c, d,  4,  7,  5)
       call transform (ff, d, a, b, c,  5, 12,  6)
       call transform (ff, c, d, a, b,  6, 17,  7)
       call transform (ff, b, c, d, a,  7, 22,  8)
       call transform (ff, a, b, c, d,  8,  7,  9)
       call transform (ff, d, a, b, c,  9, 12, 10)
       call transform (ff, c, d, a, b, 10, 17, 11)
       call transform (ff, b, c, d, a, 11, 22, 12)
       call transform (ff, a, b, c, d, 12,  7, 13)
       call transform (ff, d, a, b, c, 13, 12, 14)
       call transform (ff, c, d, a, b, 14, 17, 15)
       call transform (ff, b, c, d, a, 15, 22, 16)
       call transform (fg, a, b, c, d,  1,  5, 17)
       call transform (fg, d, a, b, c,  6,  9, 18)
       call transform (fg, c, d, a, b, 11, 14, 19)
       call transform (fg, b, c, d, a,  0, 20, 20)
       call transform (fg, a, b, c, d,  5,  5, 21)
       call transform (fg, d, a, b, c, 10,  9, 22)
       call transform (fg, c, d, a, b, 15, 14, 23)
       call transform (fg, b, c, d, a,  4, 20, 24)
       call transform (fg, a, b, c, d,  9,  5, 25)
       call transform (fg, d, a, b, c, 14,  9, 26)
       call transform (fg, c, d, a, b,  3, 14, 27)
       call transform (fg, b, c, d, a,  8, 20, 28)
       call transform (fg, a, b, c, d, 13,  5, 29)
       call transform (fg, d, a, b, c,  2,  9, 30)
       call transform (fg, c, d, a, b,  7, 14, 31)
       call transform (fg, b, c, d, a, 12, 20, 32)
       call transform (fh, a, b, c, d,  5,  4, 33)
       call transform (fh, d, a, b, c,  8, 11, 34)
       call transform (fh, c, d, a, b, 11, 16, 35)
       call transform (fh, b, c, d, a, 14, 23, 36)
       call transform (fh, a, b, c, d,  1,  4, 37)
       call transform (fh, d, a, b, c,  4, 11, 38)
       call transform (fh, c, d, a, b,  7, 16, 39)
       call transform (fh, b, c, d, a, 10, 23, 40)
       call transform (fh, a, b, c, d, 13,  4, 41)
       call transform (fh, d, a, b, c,  0, 11, 42)
       call transform (fh, c, d, a, b,  3, 16, 43)
       call transform (fh, b, c, d, a,  6, 23, 44)
       call transform (fh, a, b, c, d,  9,  4, 45)
       call transform (fh, d, a, b, c, 12, 11, 46)
       call transform (fh, c, d, a, b, 15, 16, 47)
       call transform (fh, b, c, d, a,  2, 23, 48)
       call transform (fi, a, b, c, d,  0,  6, 49)
       call transform (fi, d, a, b, c,  7, 10, 50)
       call transform (fi, c, d, a, b, 14, 15, 51)
       call transform (fi, b, c, d, a,  5, 21, 52)
       call transform (fi, a, b, c, d, 12,  6, 53)
       call transform (fi, d, a, b, c,  3, 10, 54)
       call transform (fi, c, d, a, b, 10, 15, 55)
       call transform (fi, b, c, d, a,  1, 21, 56)
       call transform (fi, a, b, c, d,  8,  6, 57)
       call transform (fi, d, a, b, c, 15, 10, 58)
       call transform (fi, c, d, a, b,  6, 15, 59)
       call transform (fi, b, c, d, a, 13, 21, 60)
       call transform (fi, a, b, c, d,  4,  6, 61)
       call transform (fi, d, a, b, c, 11, 10, 62)
       call transform (fi, c, d, a, b,  2, 15, 63)
       call transform (fi, b, c, d, a,  9, 21, 64)
       a = a + aa
       b = b + bb
       c = c + cc
       d = d + dd
       bl => bl%next
    end do
    s = digest_string ([a, b, c, d])
  contains
    subroutine transform (f, a, b, c, d, k, s, i)
      interface
         function f (x, y, z)
           import word32_t
           type(word32_t), intent(in) :: x, y, z
           type(word32_t) :: f
         end function f
      end interface
      type(word32_t), intent(inout) :: a
      type(word32_t), intent(in) :: b, c, d
      integer, intent(in) :: k, s, i
      a = b + ishftc (a + f(b, c, d) + x(k) + t(i), s)
    end subroutine transform

  end subroutine message_digest

  module function md5sum_from_string (s) result (digest)
    character(len=*), intent(in) :: s
    character(len=32) :: digest
    type(message_t) :: m
    call message_append_string (m, s)
    call message_pad (m)
    call message_digest (m, digest)
    call message_clear (m)
  end function md5sum_from_string

  module function md5sum_from_unit (u) result (digest)
    integer, intent(in) :: u
    character(len=32) :: digest
    type(message_t) :: m
    character :: char
    integer :: iostat
    READ_CHARS: do
       read (u, "(A)", advance="no", iostat=iostat)  char
       select case (iostat)
       case (0)
          call message_append_string (m, char)
       case (EOR)
          call message_append_string (m, LF)
       case (EOF)
          exit READ_CHARS
       case default
          call msg_fatal &
               ("Computing MD5 sum: I/O error while reading from scratch file")
       end select
    end do READ_CHARS
    call message_pad (m)
    call message_digest (m, digest)
    call message_clear (m)
  end function md5sum_from_unit


end submodule md5_s

