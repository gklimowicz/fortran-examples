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

submodule (cputime) cputime_s

  use io_units
  use diagnostics

  implicit none

contains

  module subroutine time_write (object, unit)
    class(time_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)", advance="no")  "Time in seconds ="
    if (object%known) then
       write (u, "(1x,ES10.3)")  object%value
    else
       write (u, "(1x,A)")  "[unknown]"
    end if
  end subroutine time_write

  module subroutine time_set_current (time)
    class(time_t), intent(out) :: time
    integer :: msecs
    call system_clock (msecs)
    time%value = real (msecs) / 1000.
    time%known = time%value > 0
  end subroutine time_set_current

  pure module subroutine real_assign_time (r, time)
    real, intent(out) :: r
    class(time_t), intent(in) :: time
    if (time%known) then
       r = time%value
    else
       r = 0
    end if
  end subroutine real_assign_time

  pure module subroutine real_default_assign_time (r, time)
    real(default), intent(out) :: r
    class(time_t), intent(in) :: time
    if (time%known) then
       r = time%value
    else
       r = 0
    end if
  end subroutine real_default_assign_time

  module subroutine time_assign_from_integer (time, ival)
    class(time_t), intent(out) :: time
    integer, intent(in) :: ival
    time%value = ival
    time%known = .true.
  end subroutine time_assign_from_integer

  module subroutine time_assign_from_real (time, rval)
    class(time_t), intent(out) :: time
    real, intent(in) :: rval
    time%value = rval
    time%known = .true.
  end subroutine time_assign_from_real

  pure module function subtract_times (t_end, t_begin) result (time)
    type(time_t) :: time
    class(time_t), intent(in) :: t_end, t_begin
    if (t_end%known .and. t_begin%known) then
       time%known = .true.
       time%value = t_end%value - t_begin%value
    end if
  end function subtract_times

  pure module function add_times (t1, t2) result (time)
    type(time_t) :: time
    class(time_t), intent(in) :: t1, t2
    if (t1%known .and. t2%known) then
       time%known = .true.
       time%value = t1%value + t2%value
    end if
  end function add_times

  module function time_is_known (time) result (flag)
    class(time_t), intent(in) :: time
    logical :: flag
    flag = time%known
  end function time_is_known

  module subroutine time_expand_s (time, sec)
    class(time_t), intent(in) :: time
    integer, intent(out) :: sec
    if (time%known) then
       sec = time%value
    else
       call msg_bug ("Time: attempt to expand undefined value")
    end if
  end subroutine time_expand_s

  module subroutine time_expand_ms (time, min, sec)
    class(time_t), intent(in) :: time
    integer, intent(out) :: min, sec
    if (time%known) then
       if (time%value >= 0) then
          sec = mod (int (time%value), 60)
       else
          sec = - mod (int (- time%value), 60)
       end if
       min = time%value / 60
    else
       call msg_bug ("Time: attempt to expand undefined value")
    end if
  end subroutine time_expand_ms

  module subroutine time_expand_hms (time, hour, min, sec)
    class(time_t), intent(in) :: time
    integer, intent(out) :: hour, min, sec
    call time%expand (min, sec)
    hour = min / 60
    if (min >= 0) then
       min = mod (min, 60)
    else
       min = - mod (-min, 60)
    end if
  end subroutine time_expand_hms

  module subroutine time_expand_dhms (time, day, hour, min, sec)
    class(time_t), intent(in) :: time
    integer, intent(out) :: day, hour, min, sec
    call time%expand (hour, min, sec)
    day = hour / 24
    if (hour >= 0) then
       hour = mod (hour, 24)
    else
       hour = - mod (- hour, 24)
    end if
  end subroutine time_expand_dhms

  module function time_to_string_s (time) result (str)
    class(time_t), intent(in) :: time
    type(string_t) :: str
    character(256) :: buffer
    integer :: s
    call time%expand (s)
    write (buffer, "(I0,'s')")  s
    str = trim (buffer)
  end function time_to_string_s

  module function time_to_string_ms (time, blank) result (str)
    class(time_t), intent(in) :: time
    logical, intent(in), optional :: blank
    type(string_t) :: str
    character(256) :: buffer
    integer :: s, m
    logical :: x_out
    x_out = .false.
    if (present (blank))  x_out = blank
    call time%expand (m, s)
    write (buffer, "(I0,'m:',I2.2,'s')")  m, abs (s)
    str = trim (buffer)
    if (x_out) then
       str = replace (str, len(str)-1, "X")
    end if
  end function time_to_string_ms

  module function time_to_string_hms (time) result (str)
    class(time_t), intent(in) :: time
    type(string_t) :: str
    character(256) :: buffer
    integer :: s, m, h
    call time%expand (h, m, s)
    write (buffer, "(I0,'h:',I2.2,'m:',I2.2,'s')")  h, abs (m), abs (s)
    str = trim (buffer)
  end function time_to_string_hms

  module function time_to_string_dhms (time) result (str)
    class(time_t), intent(in) :: time
    type(string_t) :: str
    character(256) :: buffer
    integer :: s, m, h, d
    call time%expand (d, h, m, s)
    write (buffer, "(I0,'d:',I2.2,'h:',I2.2,'m:',I2.2,'s')")  &
         d, abs (h), abs (m), abs (s)
    str = trim (buffer)
  end function time_to_string_dhms

  module subroutine timer_write (object, unit)
    class(timer_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    if (object%running) then
       write (u, "(1x,A)")  "Time in seconds = [running]"
    else
       call object%time_t%write (u)
    end if
  end subroutine timer_write

  module subroutine timer_start (timer)
    class(timer_t), intent(out) :: timer
    call timer%t1%set_current ()
    timer%running = .true.
  end subroutine timer_start

  module subroutine timer_restart (timer)
    class(timer_t), intent(inout) :: timer
    if (timer%t1%known .and. .not. timer%running) then
       timer%running = .true.
    else
       call msg_bug ("Timer: restart attempt from wrong status")
    end if
  end subroutine timer_restart

  module subroutine timer_stop (timer)
    class(timer_t), intent(inout) :: timer
    call timer%t2%set_current ()
    timer%running = .false.
    call timer%evaluate ()
  end subroutine timer_stop

  module subroutine timer_set_test_time1 (timer, t)
    class(timer_t), intent(inout) :: timer
    integer, intent(in) :: t
    timer%t1 = t
  end subroutine timer_set_test_time1

  module subroutine timer_set_test_time2 (timer, t)
    class(timer_t), intent(inout) :: timer
    integer, intent(in) :: t
    timer%t2 = t
  end subroutine timer_set_test_time2

  module subroutine timer_evaluate (timer)
    class(timer_t), intent(inout) :: timer
    timer%time_t = timer%t2 - timer%t1
  end subroutine timer_evaluate


end submodule cputime_s

