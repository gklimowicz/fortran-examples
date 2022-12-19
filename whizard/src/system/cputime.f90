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

module cputime

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: time_t
  public :: assignment(=)
  public :: timer_t

  type :: time_t
     private
     logical :: known = .false.
     real :: value = 0
   contains
     procedure :: write => time_write
     procedure :: set_current => time_set_current
     generic :: assignment(=) => time_assign_from_integer, time_assign_from_real
     procedure, private :: time_assign_from_integer
     procedure, private :: time_assign_from_real
     generic :: operator(-) => subtract_times
     generic :: operator(+) => add_times
     procedure, private :: subtract_times
     procedure, private :: add_times
     procedure :: is_known => time_is_known
     generic :: expand => time_expand_s, time_expand_ms, &
          time_expand_hms, time_expand_dhms
     procedure, private :: time_expand_s
     procedure, private :: time_expand_ms
     procedure, private :: time_expand_hms
     procedure, private :: time_expand_dhms
     procedure :: to_string_s => time_to_string_s
     procedure :: to_string_ms => time_to_string_ms
     procedure :: to_string_hms => time_to_string_hms
     procedure :: to_string_dhms => time_to_string_dhms
  end type time_t

  type, extends (time_t) :: timer_t
     private
     logical :: running = .false.
     type(time_t) :: t1, t2
   contains
     procedure :: write => timer_write
     procedure :: start => timer_start
     procedure :: restart => timer_restart
     procedure :: stop => timer_stop
     procedure :: set_test_time1 => timer_set_test_time1
     procedure :: set_test_time2 => timer_set_test_time2
     procedure :: evaluate => timer_evaluate
  end type timer_t


  interface assignment(=)
    module procedure real_assign_time
    module procedure real_default_assign_time
  end interface


  interface
    module subroutine time_write (object, unit)
      class(time_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine time_write
    module subroutine time_set_current (time)
      class(time_t), intent(out) :: time
    end subroutine time_set_current
    pure module subroutine real_assign_time (r, time)
      real, intent(out) :: r
      class(time_t), intent(in) :: time
    end subroutine real_assign_time
    pure module subroutine real_default_assign_time (r, time)
      real(default), intent(out) :: r
      class(time_t), intent(in) :: time
    end subroutine real_default_assign_time
    module subroutine time_assign_from_integer (time, ival)
      class(time_t), intent(out) :: time
      integer, intent(in) :: ival
    end subroutine time_assign_from_integer
    module subroutine time_assign_from_real (time, rval)
      class(time_t), intent(out) :: time
      real, intent(in) :: rval
    end subroutine time_assign_from_real
    pure module function subtract_times (t_end, t_begin) result (time)
      type(time_t) :: time
      class(time_t), intent(in) :: t_end, t_begin
    end function subtract_times
    pure module function add_times (t1, t2) result (time)
      type(time_t) :: time
      class(time_t), intent(in) :: t1, t2
    end function add_times
    module function time_is_known (time) result (flag)
      class(time_t), intent(in) :: time
      logical :: flag
    end function time_is_known
    module subroutine time_expand_s (time, sec)
      class(time_t), intent(in) :: time
      integer, intent(out) :: sec
    end subroutine time_expand_s
    module subroutine time_expand_ms (time, min, sec)
      class(time_t), intent(in) :: time
      integer, intent(out) :: min, sec
    end subroutine time_expand_ms
    module subroutine time_expand_hms (time, hour, min, sec)
      class(time_t), intent(in) :: time
      integer, intent(out) :: hour, min, sec
    end subroutine time_expand_hms
    module subroutine time_expand_dhms (time, day, hour, min, sec)
      class(time_t), intent(in) :: time
      integer, intent(out) :: day, hour, min, sec
    end subroutine time_expand_dhms
    module function time_to_string_s (time) result (str)
      class(time_t), intent(in) :: time
      type(string_t) :: str
    end function time_to_string_s
    module function time_to_string_ms (time, blank) result (str)
      class(time_t), intent(in) :: time
      logical, intent(in), optional :: blank
      type(string_t) :: str
    end function time_to_string_ms
    module function time_to_string_hms (time) result (str)
      class(time_t), intent(in) :: time
      type(string_t) :: str
    end function time_to_string_hms
    module function time_to_string_dhms (time) result (str)
      class(time_t), intent(in) :: time
      type(string_t) :: str
    end function time_to_string_dhms
    module subroutine timer_write (object, unit)
      class(timer_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine timer_write
    module subroutine timer_start (timer)
      class(timer_t), intent(out) :: timer
    end subroutine timer_start
    module subroutine timer_restart (timer)
      class(timer_t), intent(inout) :: timer
    end subroutine timer_restart
    module subroutine timer_stop (timer)
      class(timer_t), intent(inout) :: timer
    end subroutine timer_stop
    module subroutine timer_set_test_time1 (timer, t)
      class(timer_t), intent(inout) :: timer
      integer, intent(in) :: t
    end subroutine timer_set_test_time1
    module subroutine timer_set_test_time2 (timer, t)
      class(timer_t), intent(inout) :: timer
      integer, intent(in) :: t
    end subroutine timer_set_test_time2
    module subroutine timer_evaluate (timer)
      class(timer_t), intent(inout) :: timer
    end subroutine timer_evaluate
  end interface

end module cputime
