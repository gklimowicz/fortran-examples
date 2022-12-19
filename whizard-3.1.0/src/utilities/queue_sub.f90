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

submodule (queue) queue_s

  use, intrinsic :: iso_fortran_env, only: ERROR_UNIT

  implicit none

contains

  elemental module function queue_is_full (queue) result (flag)
    class(queue_t), intent(in) :: queue
    logical :: flag
    flag = queue%front == 1 .and. queue%rear == QUEUE_END
  end function queue_is_full

  elemental module function queue_is_empty (queue) result (flag)
    class(queue_t), intent(in) :: queue
    logical :: flag
    flag = queue%front == QUEUE_START
  end function queue_is_empty

  module subroutine queue_enqueue (queue, item)
    class(queue_t), intent(inout) :: queue
    integer, intent(in) :: item
    if (queue%is_full ()) then
       !! Do something. 
    else
       if (queue%front == QUEUE_START) queue%front = 1
       queue%rear = queue%rear + 1
       queue%item(queue%rear) = item
    end if
  end subroutine queue_enqueue

  module function queue_dequeue (queue) result (item)
    class(queue_t), intent(inout) :: queue
    integer :: item
    if (queue%is_empty ()) then
       item = 0
    else
       item = queue%item(queue%front)
       if (queue%front >= queue%rear) then
          queue%front = QUEUE_START
          queue%rear = QUEUE_START
          !! Q has only one element,
          !! so we reset the queue after deleting it.
       else
          queue%front = queue%front + 1
       end if
    end if
  end function queue_dequeue

  module function queue_peek (queue) result (item)
    class(queue_t), intent(in) :: queue
    integer :: item
    if (queue%is_empty ()) then
       item = 0
    else
       item = queue%item(queue%front)
    end if
  end function queue_peek

  module subroutine queue_write (queue, unit)
    class(queue_t), intent(in) :: queue
    integer, intent(in), optional :: unit
    integer :: u, i
    u = ERROR_UNIT; if (present (unit)) u = unit
    if (queue%is_empty ()) then
       write (u, *) "Empty Queue."
    else
       write (u, *) "Front ->", queue%front
       write (u, *) "Items ->"
       do i = 1, queue%rear
          write (u, *) queue%item(i)
       end do
       write (u, *) "Rear ->", queue%rear
    end if
  end subroutine queue_write


end submodule queue_s

