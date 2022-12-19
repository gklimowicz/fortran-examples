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

submodule (iterator) iterator_s

  use, intrinsic :: iso_fortran_env, only: ERROR_UNIT

  implicit none

contains

  module subroutine iterator_write (iter, unit)
    class(iterator_t), intent(in) :: iter
    integer, intent(in), optional :: unit
    integer :: u
    u = ERROR_UNIT; if (present (unit)) u = unit
    write (u, "(3(A,1X,I3,1X))") "CURRENT", iter%current, &
         "BEGIN", iter%begin, "END", iter%end
    flush (u)
  end subroutine iterator_write

  !! Proof: step > 0, begin < end.
  !! Proof: step < 0, begin > end.
  !! Proof: step /= 0.
  module subroutine iterator_init (iter, begin, end, step)
    class(iterator_t), intent(inout) :: iter
    integer, intent(in) :: begin
    integer, intent(in) :: end
    integer, intent(in), optional :: step
    iter%begin = begin
    iter%end = end
    iter%step = 1; if (present (step)) iter%step = step
    if (abs (iter%step) > 0) then
       iter%current = iter%begin
    else
       write (ERROR_UNIT, "(A)") "ERROR: Step size MUST be unequal to zero."
       stop 1
    end if
  end subroutine iterator_init

  pure module function iterator_at_begin (iter) result (flag)
    class(iterator_t), intent(in) :: iter
    logical :: flag
    flag = iter%current == iter%begin
  end function iterator_at_begin

  pure module function iterator_at_end (iter) result (flag)
    class(iterator_t), intent(in) :: iter
    logical :: flag
    flag = iter%current == iter%end
  end function iterator_at_end

  !! Proof: begin < current < end
  pure module function iterator_is_iterable (iter) result (flag)
    class(iterator_t), intent(in) :: iter
    logical :: flag
    if (iter%step > 0) then
       flag = iter%current <= iter%end
    else if (iter%step < 0) then
       flag = iter%current >= iter%end
    else
       flag = .false.
    end if
  end function iterator_is_iterable

  module subroutine iterator_next_step (iter)
    class(iterator_t), intent(inout) :: iter
    if (.not. iter%is_iterable ()) return
    iter%current = iter%current + iter%step
  end subroutine iterator_next_step

  !! Proof: begin <= current <= end.
  !! However, after applying the step, this does not need to be true..
  module function iterator_next (iter) result (ndx)
    class(iterator_t), intent(inout) :: iter
    integer :: ndx
    if (.not. iter%is_iterable ()) then
       ndx = 0
       return
    end if
    ndx = iter%current
    iter%current = iter%current + iter%step
  end function iterator_next

  pure module function iterator_get_current (iter) result (ndx)
    class(iterator_t), intent(in) :: iter
    integer :: ndx
    if (.not. iter%is_iterable ()) then
       ndx = 0
       return
    end if
    ndx = iter%current
  end function iterator_get_current


end submodule iterator_s

