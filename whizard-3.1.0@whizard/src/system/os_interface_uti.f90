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

module os_interface_uti

  use, intrinsic :: iso_c_binding !NODEP!

  use iso_varying_string, string_t => varying_string
  use io_units

  use os_interface

  implicit none
  private

  public :: os_interface_1

contains

  subroutine os_interface_1 (u)
    integer, intent(in) :: u
    type(dlaccess_t) :: dlaccess
    type(string_t) :: fname, libname, ext
    type(os_data_t) :: os_data
    type(string_t) :: filename_src, filename_obj
    abstract interface
       function so_test_proc (i) result (j) bind(C)
         import c_int
         integer(c_int), intent(in) :: i
         integer(c_int) :: j
       end function so_test_proc
    end interface
    procedure(so_test_proc), pointer :: so_test => null ()
    type(c_funptr) :: c_fptr
    integer :: unit
    integer(c_int) :: i
    call os_data%init ()
    fname = "so_test"
    filename_src = fname // os_data%fc_src_ext
    if (os_data%use_libtool) then
       ext = ".lo"
    else
       ext = os_data%obj_ext
    end if
    filename_obj = fname // ext
    libname = fname // '.' // os_data%fc_shrlib_ext

    write (u, "(A)")  "* Test output: OS interface"
    write (u, "(A)")  "*   Purpose: check os_interface routines"
    write (u, "(A)")

    write (u, "(A)")  "* write source file 'so_test.f90'"
    write (u, "(A)")
    unit = free_unit ()
    open (unit=unit, file=char(filename_src), action="write")
    write (unit, "(A)")  "function so_test (i) result (j) bind(C)"
    write (unit, "(A)")  "  use iso_c_binding"
    write (unit, "(A)")  "  integer(c_int), intent(in) :: i"
    write (unit, "(A)")  "  integer(c_int) :: j"
    write (unit, "(A)")  "  j = 2 * i"
    write (unit, "(A)")  "end function so_test"
    close (unit)
    write (u, "(A)")  "* compile and link as 'so_test.so/dylib'"
    write (u, "(A)")
    call os_compile_shared (fname, os_data)
    call os_link_shared (filename_obj, fname, os_data)
    write (u, "(A)")  "* load library 'so_test.so/dylib'"
    write (u, "(A)")
    call dlaccess_init (dlaccess, var_str ("."), libname, os_data)
    if (dlaccess_is_open (dlaccess)) then
       write (u, "(A)") "     success"
    else
       write (u, "(A)") "     failure"
    end if
    write (u, "(A)")  "* load symbol 'so_test'"
    write (u, "(A)")
    c_fptr = dlaccess_get_c_funptr (dlaccess, fname)
    if (c_associated (c_fptr)) then
       write (u, "(A)") "     success"
    else
       write (u, "(A)") "     failure"
    end if
    call c_f_procpointer (c_fptr, so_test)
    write (u, "(A)") "* Execute function from 'so_test.so/dylib'"
    i = 7
    write (u, "(A,1x,I1)")  "     input  = ", i
    write (u, "(A,1x,I1)")  "     result = ", so_test(i)
    if (so_test(i) / i .ne. 2) then
       write (u, "(A)")  "* Compiling and linking ISO C functions failed."
    else
       write (u, "(A)")  "* Successful."
    end if
    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"
    call dlaccess_final (dlaccess)
  end subroutine os_interface_1


end module os_interface_uti
