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
submodule (cpp_strings) cpp_strings_s

  implicit none

contains

  module subroutine cpp_string_init (s, cptr)
    class(cpp_string_t), intent(out) :: s
    type(c_ptr), intent(in) :: cptr
    s%cptr = cptr
    s%strlen = cpp_str_length (cptr)
  end subroutine cpp_string_init

  module subroutine cpp_string_final (s)
    class(cpp_string_t), intent(inout) :: s
    call cpp_str_delete (s%cptr)
    s%cptr = c_null_ptr
    s%strlen = 0
  end subroutine cpp_string_final

  module function char_from_cpp_string (s) result (c)
    type(cpp_string_t), intent(in) :: s
    character(len=s%strlen) :: c
    integer :: i
    do i = 1, s%strlen
       c(i:i) = cpp_str_get (s%cptr, int (i-1, c_int))
    end do
  end function char_from_cpp_string

  module function cpp_string_len (s) result (len)
    type(cpp_string_t), intent(in) :: s
    integer :: len
    len = s%strlen
  end function cpp_string_len

end submodule cpp_strings_s
