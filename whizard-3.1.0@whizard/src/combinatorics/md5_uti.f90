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

module md5_uti

  use diagnostics

  use md5

  implicit none
  private

  public :: md5_1

contains

  subroutine md5_1 (u)
    integer, intent(in) :: u
    character(32) :: s
    integer, parameter :: n = 7
    integer :: i
    character(80), dimension(n) :: teststring
    data teststring(1) /""/
    data teststring(2) /"a"/
    data teststring(3) /"abc"/
    data teststring(4) /"message digest"/
    data teststring(5) /"abcdefghijklmnopqrstuvwxyz"/
    data teststring(6) /"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"/
    data teststring(7) /"12345678901234567890123456789012345678901234567890123456789012345678901234567890"/
    character(32), dimension(n) :: result
    data result(1) /"D41D8CD98F00B204E9800998ECF8427E"/
    data result(2) /"0CC175B9C0F1B6A831C399E269772661"/
    data result(3) /"900150983CD24FB0D6963F7D28E17F72"/
    data result(4) /"F96B697D7CB7938D525A2F31AAF161D0"/
    data result(5) /"C3FCD3D76192E4007DFB496CCA67E13B"/
    data result(6) /"D174AB98D277D9F5A5611C2C9F419D9F"/
    data result(7) /"57EDF4A22BE3C955AC49DA2E2107B67A"/

    write (u, "(A)")  "* Test output: MD5"
    write (u, "(A)")  "*   Purpose: test MD5 sums"
    write (u, "(A)")

    do i = 1, n
       write (u, "(A)") "MD5 test string = " // '"'// &
            trim (teststring(i)) // '"'
       s = md5sum (trim (teststring(i)))
       write (u, "(A)") "MD5 check sum   = " // trim (s)
       write (u, "(A)") "Ref check sum   = " // result(i)
       if (s == result(i)) then
          call msg_message ("=> ok", u)
       else
          call msg_message ("=> MD5 sum self-test failed", u)
       end if
    end do
    call msg_message ("=============================================================================|", unit=u)
  end subroutine md5_1


end module md5_uti
