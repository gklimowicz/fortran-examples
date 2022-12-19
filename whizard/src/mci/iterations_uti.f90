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

module iterations_uti

  use iso_varying_string, string_t => varying_string

  use iterations

  implicit none
  private

  public :: iterations_1
  public :: iterations_2

contains

  subroutine iterations_1 (u)
    integer, intent(in) :: u
    type(iterations_list_t) :: it_list

    write (u, "(A)")  "* Test output: iterations_1"
    write (u, "(A)")  "*   Purpose: display empty iterations list"
    write (u, "(A)")

    call it_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: iterations_1"

  end subroutine iterations_1

  subroutine iterations_2 (u)
    integer, intent(in) :: u
    type(iterations_list_t) :: it_list

    write (u, "(A)")  "* Test output: iterations_2"
    write (u, "(A)")  "*   Purpose: fill and display iterations list"
    write (u, "(A)")

    write (u, "(A)")  "* Minimal setup (2 passes)"
    write (u, "(A)")

    call it_list%init ([2, 4], [5000, 20000])

    call it_list%write (u)
    call it_list%clear ()

    write (u, "(A)")
    write (u, "(A)")  "* Setup with flags (3 passes)"
    write (u, "(A)")

    call it_list%init ([2, 4, 5], [5000, 20000, 400], &
         [.false., .true., .true.], &
         [var_str (""), var_str ("g"), var_str ("wg")])

    call it_list%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Extract data"
    write (u, "(A)")

    write (u, "(A,I0)")  "n_pass = ", it_list%get_n_pass ()
    write (u, "(A)")
    write (u, "(A,I0)")  "n_calls(2) = ", it_list%get_n_calls (2)
    write (u, "(A)")
    write (u, "(A,I0)")  "n_it(3) = ", it_list%get_n_it (3)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: iterations_2"

  end subroutine iterations_2


end module iterations_uti
