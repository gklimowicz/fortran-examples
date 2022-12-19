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

module cputime_uti

  use iso_varying_string, string_t => varying_string

  use cputime

  implicit none
  private

  public :: cputime_1
  public :: cputime_2

contains

  subroutine cputime_1 (u)
    integer, intent(in) :: u
    type(time_t) :: time, time1, time2
    real :: t
    integer :: d, h, m, s

    write (u, "(A)")  "* Test output: cputime_1"
    write (u, "(A)")  "*   Purpose: check time operations"
    write (u, "(A)")

    write (u, "(A)") "* Undefined time"
    write (u, *)

    call time%write (u)

    write (u, *)
    write (u, "(A)") "* Set time to zero"
    write (u, *)

    time = 0
    call time%write (u)

    write (u, *)
    write (u, "(A)") "* Set time to 1.234 s"
    write (u, *)

    time = 1.234
    call time%write (u)

    t = time
    write (u, "(1x,A,F6.3)")  "Time as real =", t

    write (u, *)
    write (u, "(A)") "* Compute time difference"
    write (u, *)

    time1 = 5.33
    time2 = 7.55
    time = time2 - time1

    call time1%write (u)
    call time2%write (u)
    call time%write (u)

    write (u, *)
    write (u, "(A)") "* Compute time sum"
    write (u, *)

    time = time2 + time1

    call time1%write (u)
    call time2%write (u)
    call time%write (u)

    write (u, *)
    write (u, "(A)") "* Expand time"
    write (u, *)

    time1 = ((24 + 1) * 60 + 1) * 60 + 1
    time2 = ((3 * 24 + 23) * 60 + 59) * 60 + 59

    call time1%expand (s)
    write (u, 1)  "s =", s
    call time1%expand (m,s)
    write (u, 1)  "ms =", m, s
    call time1%expand (h,m,s)
    write (u, 1)  "hms =", h, m, s
    call time1%expand (d,h,m,s)
    write (u, 1)  "dhms =", d, h, m, s

    call time2%expand (s)
    write (u, 1)  "s =", s
    call time2%expand (m,s)
    write (u, 1)  "ms =", m, s
    call time2%expand (h,m,s)
    write (u, 1)  "hms =", h, m, s
    call time2%expand (d,h,m,s)
    write (u, 1)  "dhms =", d, h, m, s

    write (u, *)
    write (u, "(A)") "* Expand negative time"
    write (u, *)

    time1 = - (((24 + 1) * 60 + 1) * 60 + 1)
    time2 = - (((3 * 24 + 23) * 60 + 59) * 60 + 59)

    call time1%expand (s)
    write (u, 1)  "s =", s
    call time1%expand (m,s)
    write (u, 1)  "ms =", m, s
    call time1%expand (h,m,s)
    write (u, 1)  "hms =", h, m, s
    call time1%expand (d,h,m,s)
    write (u, 1)  "dhms =", d, h, m, s

    call time2%expand (s)
    write (u, 1)  "s =", s
    call time2%expand (m,s)
    write (u, 1)  "ms =", m, s
    call time2%expand (h,m,s)
    write (u, 1)  "hms =", h, m, s
    call time2%expand (d,h,m,s)
    write (u, 1)  "dhms =", d, h, m, s

1   format (1x,A,1x,4(I0,:,':'))

    write (u, *)
    write (u, "(A)") "* String from time"
    write (u, *)

    time1 = ((24 + 1) * 60 + 1) * 60 + 1
    time2 = ((3 * 24 + 23) * 60 + 59) * 60 + 59

    write (u, "(A)")  char (time1%to_string_s ())
    write (u, "(A)")  char (time1%to_string_ms ())
    write (u, "(A)")  char (time1%to_string_hms ())
    write (u, "(A)")  char (time1%to_string_dhms ())

    write (u, "(A)")  char (time2%to_string_s ())
    write (u, "(A)")  char (time2%to_string_ms ())
    write (u, "(A)")  char (time2%to_string_hms ())
    write (u, "(A)")  char (time2%to_string_dhms ())

    write (u, "(A)")
    write (u, "(A)")  "* Blanking out the last second entry"
    write (u, "(A)")

    write (u, "(A)")  char (time1%to_string_ms ())
    write (u, "(A)")  char (time1%to_string_ms (.true.))

    write (u, *)
    write (u, "(A)") "* String from negative time"
    write (u, *)

    time1 = -(((24 + 1) * 60 + 1) * 60 + 1)
    time2 = -(((3 * 24 + 23) * 60 + 59) * 60 + 59)

    write (u, "(A)")  char (time1%to_string_s ())
    write (u, "(A)")  char (time1%to_string_ms ())
    write (u, "(A)")  char (time1%to_string_hms ())
    write (u, "(A)")  char (time1%to_string_dhms ())

    write (u, "(A)")  char (time2%to_string_s ())
    write (u, "(A)")  char (time2%to_string_ms ())
    write (u, "(A)")  char (time2%to_string_hms ())
    write (u, "(A)")  char (time2%to_string_dhms ())

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: cputime_1"

  end subroutine cputime_1

  subroutine cputime_2 (u)
    integer, intent(in) :: u
    type(timer_t) :: timer

    write (u, "(A)")  "* Test output: cputime_2"
    write (u, "(A)")  "*   Purpose: check timer"
    write (u, "(A)")

    write (u, "(A)") "* Undefined timer"
    write (u, *)

    call timer%write (u)

    write (u, *)
    write (u, "(A)") "* Start timer"
    write (u, *)

    call timer%start ()
    call timer%write (u)

    write (u, *)
    write (u, "(A)") "* Stop timer (injecting fake timings)"
    write (u, *)

    call timer%stop ()
    call timer%set_test_time1 (2)
    call timer%set_test_time2 (5)
    call timer%evaluate ()
    call timer%write (u)

    write (u, *)
    write (u, "(A)") "* Restart timer"
    write (u, *)

    call timer%restart ()
    call timer%write (u)

    write (u, *)
    write (u, "(A)") "* Stop timer again (injecting fake timing)"
    write (u, *)

    call timer%stop ()
    call timer%set_test_time2 (10)
    call timer%evaluate ()
    call timer%write (u)

    write (u, *)
    write (u, "(A)")  "* Test output end: cputime_2"

  end subroutine cputime_2


end module cputime_uti
