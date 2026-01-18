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

module cascades2_lexer_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use numeric_utils

  use cascades2_lexer

  implicit none
  private

  public :: cascades2_lexer_1

contains

  subroutine cascades2_lexer_1 (u)
    integer, intent(in) :: u
    integer :: u_in = 8
    character(len=300) :: line
    integer :: stat
    logical :: fail
    type(dag_string_t) :: dag_string

    write (u, "(A)")  "* Test output: cascades2_lexer_1"
    write (u, "(A)")  "*   Purpose: read lines of O'Mega's phase space output, translate"
    write (u, "(A)")  "*            to dag_string, retranslate to character string and"
    write (u, "(A)")  "*            compare"
    write (u, "(A)")

    open (unit=u_in, file="cascades2_lexer_1.fds", status='old', action='read')

    stat = 0
    fail = .false.
    read (unit=u_in, fmt="(A)", iostat=stat) line
    do while (stat == 0 .and. .not. fail)
       read (unit=u_in, fmt="(A)", iostat=stat) line
       if (stat /= 0) exit
       dag_string = line
       fail = (char(dag_string) /= line)
    end do
    if (fail) then
       write (u, "(A)")  "* Test result: Test failed!"
    else
       write (u, "(A)")  "* Test result: Test passed"
    end if

    close (u_in)
    write (u, *)
    write (u, "(A)")  "* Test output end: cascades2_lexer_1"
  end subroutine cascades2_lexer_1


end module cascades2_lexer_uti
