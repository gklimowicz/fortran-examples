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

module queue

  implicit none
  private

  public :: queue_t

  integer, parameter :: QUEUE_SIZE = 10, &
       QUEUE_START = 0, &
       QUEUE_END = QUEUE_SIZE


  type :: queue_t
     private
     integer, dimension(QUEUE_SIZE) :: item
     integer :: front = 0
     integer :: rear = 0
   contains
     procedure :: is_full => queue_is_full
     procedure :: is_empty => queue_is_empty
     procedure :: enqueue => queue_enqueue
     procedure :: dequeue => queue_dequeue
     procedure :: peek => queue_peek
     procedure :: write => queue_write
  end type queue_t


  interface
    elemental module function queue_is_full (queue) result (flag)
      class(queue_t), intent(in) :: queue
      logical :: flag
    end function queue_is_full
    elemental module function queue_is_empty (queue) result (flag)
      class(queue_t), intent(in) :: queue
      logical :: flag
    end function queue_is_empty
    module subroutine queue_enqueue (queue, item)
      class(queue_t), intent(inout) :: queue
      integer, intent(in) :: item
    end subroutine queue_enqueue
    module function queue_dequeue (queue) result (item)
      class(queue_t), intent(inout) :: queue
      integer :: item
    end function queue_dequeue
    module function queue_peek (queue) result (item)
      class(queue_t), intent(in) :: queue
      integer :: item
    end function queue_peek
    module subroutine queue_write (queue, unit)
      class(queue_t), intent(in) :: queue
      integer, intent(in), optional :: unit
    end subroutine queue_write
  end interface

end module queue
