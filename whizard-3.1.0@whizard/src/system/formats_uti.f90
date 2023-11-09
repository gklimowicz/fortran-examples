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

module formats_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  use formats

  implicit none
  private

  public :: format_1



contains

  subroutine format_1 (u)
    integer, intent(in) :: u
    write (u, "(A)")  "*** Test 1: a string ***"
    write (u, "(A)")
    call test_run (var_str("%s"), 1, [4], ['abcdefghij'], u)
    write (u, "(A)")  "*** Test 2: two integers ***"
    write (u, "(A)")
    call test_run (var_str("%d,%d"), 2, [2, 2], ['42', '13'], u)
    write (u, "(A)")  "*** Test 3: floating point number ***"
    write (u, "(A)")
    call test_run (var_str("%8.4f"), 1, [3], ['42567.12345'], u)
    write (u, "(A)")  "*** Test 4: general expression ***"
    call test_run (var_str("%g"), 1, [3], ['3.1415'], u)
    contains
      subroutine test_run (fmt, n_args, type, buffer, unit)
        type(string_t), intent(in) :: fmt
        integer, intent(in) :: n_args, unit
        logical :: lval
        integer :: ival
        real(default) :: rval
        integer :: i
        type(string_t) :: string
        type(sprintf_arg_t), dimension(:), allocatable :: arg
        integer, dimension(n_args), intent(in) :: type
        character(*), dimension(n_args), intent(in) :: buffer
        write (unit, "(A,A)")   "Format string :", char(fmt)
        write (unit, "(A,I1)")  "Number of args:", n_args
        allocate (arg (n_args))
        do i = 1, n_args
           write (unit, "(A,I1)")  "Argument (type ) = ", type(i)
           select case (type(i))
           case (ARGTYPE_LOG)
              read (buffer(i), *)  lval
              call sprintf_arg_init (arg(i), lval)
           case (ARGTYPE_INT)
              read (buffer(i), *)  ival
              call sprintf_arg_init (arg(i), ival)
           case (ARGTYPE_REAL)
              read (buffer(i), *)  rval
              call sprintf_arg_init (arg(i), rval)
           case (ARGTYPE_STR)
              call sprintf_arg_init (arg(i), var_str (trim (buffer(i))))
           end select
         end do
         string = sprintf (fmt, arg)
         write (unit, "(A,A,A)")  "Result: '", char (string), "'"
         deallocate (arg)
       end subroutine test_run
  end subroutine format_1


end module formats_uti
