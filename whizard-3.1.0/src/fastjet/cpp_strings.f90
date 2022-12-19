!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

! Wrapper module for C++ string handling
! It defines a type cpp_string_t that acts as a proxy to a C++ string
module cpp_strings

  use, intrinsic :: iso_c_binding

  implicit none
  private

  public :: cpp_string_t
  public :: char, len

  type :: cpp_string_t
     private
     type(c_ptr) :: cptr = c_null_ptr
     integer :: strlen = 0
   contains
     procedure :: init => cpp_string_init
     procedure :: final => cpp_string_final
  end type cpp_string_t

  interface
     subroutine cpp_str_delete (cpp_str) bind (C)
       import
       type(c_ptr), value :: cpp_str
     end subroutine cpp_str_delete
  end interface

  interface
     function cpp_str_length (cpp_str) bind (C) result (length)
       import
       type(c_ptr), intent(in), value :: cpp_str
       integer(c_int) :: length
     end function cpp_str_length
  end interface
       
  interface
     function cpp_str_get (cpp_str, i) bind (C) result (c)
       import
       type(c_ptr), intent(in), value :: cpp_str
       integer(c_int), intent(in), value :: i
       character(c_char) :: c
     end function cpp_str_get
  end interface
       
  interface char
     module procedure char_from_cpp_string
  end interface char

  interface len
     module procedure cpp_string_len
  end interface len

  interface
    module subroutine cpp_string_init (s, cptr)
      class(cpp_string_t), intent(out) :: s
      type(c_ptr), intent(in) :: cptr
    end subroutine cpp_string_init
    module subroutine cpp_string_final (s)
      class(cpp_string_t), intent(inout) :: s
    end subroutine cpp_string_final
    module function char_from_cpp_string (s) result (c)
      type(cpp_string_t), intent(in) :: s
      character(len=s%strlen) :: c
    end function char_from_cpp_string
    module function cpp_string_len (s) result (len)
      type(cpp_string_t), intent(in) :: s
      integer :: len
    end function cpp_string_len
  end interface

end module cpp_strings
