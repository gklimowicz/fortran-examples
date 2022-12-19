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

submodule (formats) formats_s

  use io_units
  use diagnostics

  implicit none

contains

  module subroutine sprintf_arg_init_log (arg, lval)
    type(sprintf_arg_t), intent(out) :: arg
    logical, intent(in) :: lval
    arg%type = ARGTYPE_STR
    if (lval) then
       allocate (arg%sval (5))
       arg%sval = ['t', 'r', 'u', 'e', c_null_char]
    else
       allocate (arg%sval (6))
       arg%sval = ['f', 'a', 'l', 's', 'e', c_null_char]
    end if
  end subroutine sprintf_arg_init_log

  module subroutine sprintf_arg_init_int (arg, ival)
    type(sprintf_arg_t), intent(out) :: arg
    integer, intent(in) :: ival
    arg%type = ARGTYPE_INT
    allocate (arg%ival (1))
    arg%ival = ival
  end subroutine sprintf_arg_init_int

  module subroutine sprintf_arg_init_real (arg, rval)
    type(sprintf_arg_t), intent(out) :: arg
    real(default), intent(in) :: rval
    arg%type = ARGTYPE_REAL
    allocate (arg%rval (1))
    arg%rval = rval
  end subroutine sprintf_arg_init_real

  module subroutine sprintf_arg_init_str (arg, sval)
    type(sprintf_arg_t), intent(out) :: arg
    type(string_t), intent(in) :: sval
    integer :: i
    arg%type = ARGTYPE_STR
    allocate (arg%sval (len (sval) + 1))
    do i = 1, len (sval)
       arg%sval(i) = extract (sval, i, i)
    end do
    arg%sval(len (sval) + 1) = c_null_char
  end subroutine sprintf_arg_init_str

  subroutine sprintf_arg_write (arg, unit)
    type(sprintf_arg_t), intent(in) :: arg
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    select case (arg%type)
    case (ARGTYPE_NONE)
      write (u, *) "[none]"
    case (ARGTYPE_INT)
      write (u, "(1x,A,1x)", advance = "no")  "[int]"
      write (u, *)  arg%ival
    case (ARGTYPE_REAL)
      write (u, "(1x,A,1x)", advance = "no")  "[real]"
      write (u, *)  arg%rval
    case (ARGTYPE_STR)
      write (u, "(1x,A,1x,A)", advance = "no")  "[string]", '"'
      write (u, *)  arg%rval, '"'
    end select
  end subroutine sprintf_arg_write

  elemental function sprintf_arg_get_length (arg) result (length)
    integer :: length
    type(sprintf_arg_t), intent(in) :: arg
    select case (arg%type)
    case (ARGTYPE_INT)
       length = log10 (real (huge (arg%ival(1)))) + 2
    case (ARGTYPE_REAL)
       length = log10 (real (radix (arg%rval(1))) ** digits (arg%rval(1))) + 8
    case (ARGTYPE_STR)
       length = size (arg%sval)
    case default
       length = 0
    end select
  end function sprintf_arg_get_length

  subroutine sprintf_arg_apply_sprintf (arg, fmt, result, actual_length)
    type(sprintf_arg_t), intent(in) :: arg
    character(c_char), dimension(:), intent(in) :: fmt
    character(c_char), dimension(:), intent(inout) :: result
    integer, intent(out) :: actual_length
    integer(c_int) :: ival
    real(c_double) :: rval
    select case (arg%type)
    case (ARGTYPE_NONE)
      actual_length = sprintf_none (result, fmt)
    case (ARGTYPE_INT)
      ival = arg%ival(1)
      actual_length = sprintf_int (result, fmt, ival)
    case (ARGTYPE_REAL)
      rval = arg%rval(1)
      actual_length = sprintf_double (result, fmt, rval)
    case (ARGTYPE_STR)
      actual_length = sprintf_str (result, fmt, arg%sval)
    case default
      call msg_bug ("sprintf_arg_apply_sprintf called with illegal type")
    end select
    if (actual_length < 0) then
       write (msg_buffer, *) "Format: '", fmt, "'"
       call msg_message ()
       write (msg_buffer, *) "Output: '", result, "'"
       call msg_message ()
       call msg_error ("I/O error in sprintf call")
       actual_length = 0
    end if
  end subroutine sprintf_arg_apply_sprintf

  subroutine sprintf_interface_init (intf, fmt, arg)
    type(sprintf_interface_t), intent(out) :: intf
    type(string_t), intent(in) :: fmt
    type(sprintf_arg_t), intent(in) :: arg
    integer :: fmt_len, i
    fmt_len = len (fmt)
    allocate (intf%input_fmt (fmt_len + 1))
    do i = 1, fmt_len
       intf%input_fmt(i) = extract (fmt, i, i)
    end do
    intf%input_fmt(fmt_len+1) = c_null_char
    intf%arg = arg
    allocate (intf%output_str (len (fmt) + sprintf_arg_get_length (arg) + 1))
  end subroutine sprintf_interface_init

  subroutine sprintf_interface_write (intf, unit)
    type(sprintf_interface_t), intent(in) :: intf
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, *) "Format string = ", '"', intf%input_fmt, '"'
    write (u, "(1x,A,1x)", advance = "no") "Argument = "
    call sprintf_arg_write (intf%arg, unit)
    if (intf%output_str_len > 0) then
       write (u, *) "Result string = ", &
            '"', intf%output_str (1:intf%output_str_len), '"'
    end if
  end subroutine sprintf_interface_write

  function sprintf_interface_get_result (intf) result (string)
    type(string_t) :: string
    type(sprintf_interface_t), intent(in) :: intf
    character(kind = c_char, len = max (intf%output_str_len, 0)) :: buffer
    integer :: i
    if (intf%output_str_len > 0) then
       do i = 1, intf%output_str_len
          buffer(i:i) = intf%output_str(i)
       end do
       string = buffer(1:intf%output_str_len)
    else
       string = ""
    end if
  end function sprintf_interface_get_result

  subroutine sprintf_interface_apply_sprintf (intf)
    type(sprintf_interface_t), intent(inout) :: intf
    call sprintf_arg_apply_sprintf &
         (intf%arg, intf%input_fmt, intf%output_str, intf%output_str_len)
  end subroutine sprintf_interface_apply_sprintf

  subroutine chop_and_check_format_string (fmt, arg, intf)
    type(string_t), intent(in) :: fmt
    type(sprintf_arg_t), dimension(:), intent(in) :: arg
    type(sprintf_interface_t), dimension(:), intent(out), allocatable :: intf
    integer :: n_args, i
    type(string_t), dimension(:), allocatable :: split_fmt
    type(string_t) :: word, buffer, separator
    integer :: pos, length, l
    logical :: ok
    type(sprintf_arg_t) :: arg_null
    ok = .true.
    length = 0
    n_args = size (arg)
    allocate (split_fmt (0:n_args))
    split_fmt = ""
    buffer = fmt
    SCAN_ARGS: do i = 1, n_args
       FIND_CONVERSION: do
          call split (buffer, word, "%", separator=separator)
          if (separator == "") then
             call msg_message ('"' // char (fmt) // '"')
             call msg_error ("C-formatting string: " &
                  // "too few conversion specifiers in format string")
             ok = .false.;  exit SCAN_ARGS
          end if
          split_fmt(i-1) = split_fmt(i-1) // word
          if (extract (buffer, 1, 1) /= "%") then
             split_fmt(i) = "%"
             exit FIND_CONVERSION
          else
             split_fmt(i-1) = split_fmt(i-1) // "%"
          end if
       end do FIND_CONVERSION
       pos = verify (buffer, "#0-+ ")   ! Flag characters (zero or more)
       split_fmt(i) = split_fmt(i) // extract (buffer, 1, pos-1)
       buffer = remove (buffer, 1, pos-1)
       pos = verify (buffer, "123456890")  ! Field width
       word = extract (buffer, 1, pos-1)
       if (len (word) /= 0) then
         call read_int_from_string (word, len (word), l)
         length = length + l
       end if
       split_fmt(i) = split_fmt(i) // word
       buffer = remove (buffer, 1, pos-1)
       if (extract (buffer, 1, 1) == ".") then
          buffer = remove (buffer, 1, 1)
          pos = verify (buffer, "1234567890")   ! Precision
          split_fmt(i) = split_fmt(i) // "." // extract (buffer, 1, pos-1)
          buffer = remove (buffer, 1, pos-1)
       end if
       ! Length modifier would come here, but is not allowed
       select case (char (extract (buffer, 1, 1)))  ! conversion specifier
       case ("d", "i")
          if (arg(i)%type /= ARGTYPE_INT) then
             call msg_message ('"' // char (fmt) // '"')
             call msg_error ("C-formatting string: " &
                  // "argument type mismatch: integer value expected")
             ok = .false.;  exit SCAN_ARGS
          end if
       case ("e", "E", "f", "F", "g", "G")
          if (arg(i)%type /= ARGTYPE_REAL) then
             call msg_message ('"' // char (fmt) // '"')
             call msg_error ("C-formatting string: " &
                  // "argument type mismatch: real value expected")
             ok = .false.;  exit SCAN_ARGS
          end if
       case ("s")
          if (arg(i)%type /= ARGTYPE_STR) then
             call msg_message ('"' // char (fmt) // '"')
             call msg_error ("C-formatting string: " &
                  // "argument type mismatch: logical or string value expected")
             ok = .false.;  exit SCAN_ARGS
          end if
       case default
          call msg_message ('"' // char (fmt) // '"')
          call msg_error ("C-formatting string: " &
               // "illegal or incomprehensible conversion specifier")
          ok = .false.;  exit SCAN_ARGS
       end select
       split_fmt(i) = split_fmt(i) // extract (buffer, 1, 1)
       buffer = remove (buffer, 1, 1)
    end do SCAN_ARGS
    if (ok) then
       FIND_EXTRA_CONVERSION: do
          call split (buffer, word, "%", separator=separator)
          split_fmt(n_args) = split_fmt(n_args) // word // separator
          if (separator == "")  exit FIND_EXTRA_CONVERSION
          if (extract (buffer, 1, 1) == "%") then
             split_fmt(n_args) = split_fmt(n_args) // "%"
             buffer = remove (buffer, 1, 1)
          else
             call msg_message ('"' // char (fmt) // '"')
             call msg_error ("C-formatting string: " &
                  // "too many conversion specifiers in format string")
             ok = .false.;  exit FIND_EXTRA_CONVERSION
          end if
       end do FIND_EXTRA_CONVERSION
       split_fmt(n_args) = split_fmt(n_args) // buffer
       allocate (intf (0:n_args))
       call sprintf_interface_init (intf(0), split_fmt(0), arg_null)
       do i = 1, n_args
          call sprintf_interface_init (intf(i), split_fmt(i), arg(i))
       end do
    else
       allocate (intf (0))
    end if
  contains
    subroutine read_int_from_string (word, length, l)
      type(string_t), intent(in) :: word
      integer, intent(in) :: length
      integer, intent(out) :: l
      character(len=length) :: buffer
      buffer = word
      read (buffer, *) l
    end subroutine read_int_from_string
  end subroutine chop_and_check_format_string

  module function sprintf (fmt, arg) result (string)
    type(string_t) :: string
    type(string_t), intent(in) :: fmt
    type(sprintf_arg_t), dimension(:), intent(in) :: arg
    type(sprintf_interface_t), dimension(:), allocatable :: intf
    integer :: i
    string = ""
    call chop_and_check_format_string (fmt, arg, intf)
    if (size (intf) > 0) then
       do i = 0, ubound (intf, 1)
          call sprintf_interface_apply_sprintf (intf(i))
          string = string // sprintf_interface_get_result (intf(i))
       end do
    end if
  end function sprintf


end submodule formats_s

