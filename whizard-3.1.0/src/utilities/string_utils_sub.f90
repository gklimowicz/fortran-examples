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

submodule (string_utils) string_utils_s

  implicit none

contains

  module function upper_case_char (string) result (new_string)
    character(*), intent(in) :: string
    character(len(string)) :: new_string
    integer :: pos, code
    integer, parameter :: offset = ichar('A')-ichar('a')
    do pos = 1, len (string)
       code = ichar (string(pos:pos))
       select case (code)
       case (ichar('a'):ichar('z'))
          new_string(pos:pos) = char (code + offset)
       case default
          new_string(pos:pos) = string(pos:pos)
       end select
    end do
  end function upper_case_char

  module function lower_case_char (string) result (new_string)
    character(*), intent(in) :: string
    character(len(string)) :: new_string
    integer :: pos, code
    integer, parameter :: offset = ichar('a')-ichar('A')
    do pos = 1, len (string)
       code = ichar (string(pos:pos))
       select case (code)
       case (ichar('A'):ichar('Z'))
          new_string(pos:pos) = char (code + offset)
       case default
          new_string(pos:pos) = string(pos:pos)
       end select
    end do
  end function lower_case_char

  module function upper_case_string (string) result (new_string)
    type(string_t), intent(in) :: string
    type(string_t) :: new_string
    new_string = upper_case_char (char (string))
  end function upper_case_string

  module function lower_case_string (string) result (new_string)
    type(string_t), intent(in) :: string
    type(string_t) :: new_string
    new_string = lower_case_char (char (string))
  end function lower_case_string

  pure module function string_f2c_char (i) result (o)
    character(*), intent(in) :: i
    character(kind=c_char, len=len (i) + 1) :: o
    o = i // c_null_char
  end function string_f2c_char

  pure module function string_f2c_var_str (i) result (o)
    type(string_t), intent(in) :: i
    character(kind=c_char, len=len (i) + 1) :: o
    o = char (i) // c_null_char
  end function string_f2c_var_str

  module subroutine strcpy_f2c (fstring, cstring)
    character(*), intent(in) :: fstring
    character(c_char), dimension(*), intent(inout) :: cstring
    
    integer :: i

    do i = 1, len (fstring)
       cstring(i) = fstring(i:i)
    end do
    cstring(len(fstring)+1) = c_null_char
    
  end subroutine strcpy_f2c

  module function string_c2f (cstring) result (fstring)
    character(c_char), dimension(*), intent(in) :: cstring
    character(:), allocatable :: fstring
    
    integer :: i, n
    
    n = 0
    do while (cstring(n+1) /= c_null_char)
       n = n + 1
    end do
    
    allocate (character(n) :: fstring)
    do i = 1, n
       fstring(i:i) = cstring(i)
    end do

  end function string_c2f

  module function str_log (l) result (s)
    logical, intent(in) :: l
    type(string_t) :: s
    if (l) then
       s = "True"
    else
       s = "False"
    end if
  end function str_log

  module function str_logs (x) result (s)
    logical, dimension(:), intent(in) :: x
    type(string_t) :: s
    integer :: i
    s = '['
    do i = 1, size(x) - 1
       s = s // str(x(i)) // ', '
    end do
    s = s // str(x(size(x))) // ']'
  end function str_logs

  module function str_int (i) result (s)
    integer, intent(in) :: i
    type(string_t) :: s
    character(32) :: buffer
    write (buffer, "(I0)")  i
    s = var_str (trim (adjustl (buffer)))
  end function str_int

  module function str_ints (x) result (s)
    integer, dimension(:), intent(in) :: x
    type(string_t) :: s
    integer :: i
    s = '['
    do i = 1, size(x) - 1
       s = s // str(x(i)) // ', '
    end do
    s = s // str(x(size(x))) // ']'
  end function str_ints

  module function str_real (x) result (s)
    real(default), intent(in) :: x
    type(string_t) :: s
    character(32) :: buffer
    write (buffer, "(ES17.10)")  x
    s = var_str (trim (adjustl (buffer)))
  end function str_real

  module function str_reals (x) result (s)
    real(default), dimension(:), intent(in) :: x
    type(string_t) :: s
    integer :: i
    s = '['
    do i = 1, size(x) - 1
       s = s // str(x(i)) // ', '
    end do
    s = s // str(x(size(x))) // ']'
  end function str_reals

  module function str_complex (x) result (s)
    complex(default), intent(in) :: x
    type(string_t) :: s
    s = str_real (real (x)) // " + i " // str_real (aimag (x))
  end function str_complex

  module function str_complexs (x) result (s)
    complex(default), dimension(:), intent(in) :: x
    type(string_t) :: s
    integer :: i
    s = '['
    do i = 1, size(x) - 1
       s = s // str(x(i)) // ', '
    end do
    s = s // str(x(size(x))) // ']'
  end function str_complexs

  module function read_rval (s) result (rval)
    real(default) :: rval
    type(string_t), intent(in) :: s
    character(80) :: buffer
    buffer = s
    read (buffer, *)  rval
  end function read_rval

  module function read_ival (s) result (ival)
    integer :: ival
    type(string_t), intent(in) :: s
    character(80) :: buffer
    buffer = s
    read (buffer, *)  ival
  end function read_ival

  pure module function string_contains_word &
       (str, word, include_identical) result (val)
    logical :: val
    type(string_t), intent(in) :: str, word
    type(string_t) :: str_tmp, str_out
    logical, intent(in), optional :: include_identical
    logical :: yorn
    str_tmp = str
    val = .false.
    yorn = .false.; if (present (include_identical))  yorn = include_identical
    if (yorn)  val = str == word
    call split (str_tmp, str_out, word)
    val = val .or. (str_out /= "")
  end function string_contains_word

  pure module subroutine split_string (str, separator, str_array)
    type(string_t), dimension(:), allocatable, intent(out) :: str_array
    type(string_t), intent(in) :: str, separator
    type(string_t) :: str_tmp, str_out
    integer :: n_str
    n_str = 0; str_tmp = str
    do while (string_contains_word (str_tmp, separator))
       n_str = n_str + 1
       call split (str_tmp, str_out, separator)
    end do
    allocate (str_array (n_str))
    n_str = 1; str_tmp = str
    do while (string_contains_word (str_tmp, separator))
       call split (str_tmp, str_array (n_str), separator)
       n_str = n_str + 1
    end do
  end subroutine split_string


end submodule string_utils_s

