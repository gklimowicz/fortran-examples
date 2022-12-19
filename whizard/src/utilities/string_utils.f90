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

module string_utils

  use, intrinsic :: iso_c_binding

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: upper_case
  public :: lower_case
  public :: string_f2c
  public :: strcpy_f2c
  public :: string_c2f
  public :: str
  public :: read_rval
  public :: read_ival
  public :: string_contains_word
  public :: split_string

  interface upper_case
     module procedure upper_case_char, upper_case_string
  end interface
  interface lower_case
     module procedure lower_case_char, lower_case_string
  end interface
  interface string_f2c
     module procedure string_f2c_char, string_f2c_var_str
  end interface string_f2c
  interface str
     module procedure str_log, str_logs, str_int, str_ints, &
            str_real, str_reals, str_complex, str_complexs
  end interface

  interface
    module function upper_case_char (string) result (new_string)
      character(*), intent(in) :: string
      character(len(string)) :: new_string
    end function upper_case_char
    module function lower_case_char (string) result (new_string)
      character(*), intent(in) :: string
      character(len(string)) :: new_string
    end function lower_case_char
    module function upper_case_string (string) result (new_string)
      type(string_t), intent(in) :: string
      type(string_t) :: new_string
    end function upper_case_string
    module function lower_case_string (string) result (new_string)
      type(string_t), intent(in) :: string
      type(string_t) :: new_string
    end function lower_case_string
    pure module function string_f2c_char (i) result (o)
      character(*), intent(in) :: i
      character(kind=c_char, len=len (i) + 1) :: o
    end function string_f2c_char
    pure module function string_f2c_var_str (i) result (o)
      type(string_t), intent(in) :: i
      character(kind=c_char, len=len (i) + 1) :: o
    end function string_f2c_var_str
    module subroutine strcpy_f2c (fstring, cstring)
      character(*), intent(in) :: fstring
      character(c_char), dimension(*), intent(inout) :: cstring
    end subroutine strcpy_f2c
    module function string_c2f (cstring) result (fstring)
      character(c_char), dimension(*), intent(in) :: cstring
      character(:), allocatable :: fstring
    end function string_c2f
    module function str_log (l) result (s)
      logical, intent(in) :: l
      type(string_t) :: s
    end function str_log
    module function str_logs (x) result (s)
      logical, dimension(:), intent(in) :: x
      type(string_t) :: s
    end function str_logs
    module function str_int (i) result (s)
      integer, intent(in) :: i
      type(string_t) :: s
    end function str_int
    module function str_ints (x) result (s)
      integer, dimension(:), intent(in) :: x
      type(string_t) :: s
    end function str_ints
    module function str_real (x) result (s)
      real(default), intent(in) :: x
      type(string_t) :: s
    end function str_real
    module function str_reals (x) result (s)
      real(default), dimension(:), intent(in) :: x
      type(string_t) :: s
    end function str_reals
    module function str_complex (x) result (s)
      complex(default), intent(in) :: x
      type(string_t) :: s
    end function str_complex
    module function str_complexs (x) result (s)
      complex(default), dimension(:), intent(in) :: x
      type(string_t) :: s
    end function str_complexs
    module function read_rval (s) result (rval)
      real(default) :: rval
      type(string_t), intent(in) :: s
    end function read_rval
    module function read_ival (s) result (ival)
      integer :: ival
      type(string_t), intent(in) :: s
    end function read_ival
    pure module function string_contains_word &
         (str, word, include_identical) result (val)
      logical :: val
      type(string_t), intent(in) :: str, word
      logical, intent(in), optional :: include_identical
    end function string_contains_word
    pure module subroutine split_string (str, separator, str_array)
      type(string_t), dimension(:), allocatable, intent(out) :: str_array
      type(string_t), intent(in) :: str, separator
    end subroutine split_string
  end interface

end module string_utils
