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
module iterator

  implicit none
  private

  public :: iterator_t

  !! Forward
  type :: iterator_t
     integer :: current = 0
     integer :: begin = 0
     integer :: end = 0
     integer :: step = 1
   contains
     procedure :: write => iterator_write
     procedure :: init => iterator_init
     procedure :: at_begin => iterator_at_begin
     procedure :: at_end => iterator_at_end
     procedure :: is_iterable => iterator_is_iterable
     procedure :: next_step => iterator_next_step
     procedure :: next => iterator_next
     procedure :: get_current => iterator_get_current
  end type iterator_t


  interface
    module subroutine iterator_write (iter, unit)
      class(iterator_t), intent(in) :: iter
      integer, intent(in), optional :: unit
    end subroutine iterator_write
   module subroutine iterator_init (iter, begin, end, step)
     class(iterator_t), intent(inout) :: iter
     integer, intent(in) :: begin
     integer, intent(in) :: end
     integer, intent(in), optional :: step
   end subroutine iterator_init
    pure module function iterator_at_begin (iter) result (flag)
      class(iterator_t), intent(in) :: iter
      logical :: flag
    end function iterator_at_begin
    pure module function iterator_at_end (iter) result (flag)
      class(iterator_t), intent(in) :: iter
      logical :: flag
    end function iterator_at_end
    pure module function iterator_is_iterable (iter) result (flag)
      class(iterator_t), intent(in) :: iter
      logical :: flag
    end function iterator_is_iterable
    module subroutine iterator_next_step (iter)
      class(iterator_t), intent(inout) :: iter
    end subroutine iterator_next_step
    module function iterator_next (iter) result (ndx)
      class(iterator_t), intent(inout) :: iter
      integer :: ndx
    end function iterator_next
    pure module function iterator_get_current (iter) result (ndx)
      class(iterator_t), intent(in) :: iter
      integer :: ndx
    end function iterator_get_current
  end interface

end module iterator

