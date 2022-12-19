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

module array_list_uti

  use array_list

  implicit none
  private

public :: array_list_1

 contains

subroutine array_list_1 (u)
  integer, intent(in) :: u
  type(array_list_t) :: list
  integer :: ndx, data

  write (u, "(A)") "* Test output: Array list"
  write (u, "(A)") "*   Purpose: test interface and implementation of array list"
  write (u, "(A)")

  write (u, "(A)") "* Init array_list_t ..."
  call list%init ()
  write (u, "(A)") "* Test adding a single element..."
  call list%add (1)
  write (u, "(A)") "* Test removing a single element..."
  data = list%remove ()
  write (u, "(A)") "* Test growing (unnecessary, so just return)..."
  call list%grow_size ()
  write (u, "(A)") "* Test adding elements beyond initial capacity..."
  call test_grow_and_add (list)
  write (u, "(A)") "* Test adding at specific position..."
  call list%add_at (10, -1)
  write (u, "(A)") "* Test removing at specific position..."
  data = list%remove_at (11)
  write (u, "(A)") "* Test reverse ordering..."
  call list%reverse_order ()
  write (u, "(A)") "* Test sorting..."
  call list%sort ()
  write (u, "(A)") "* Test finding..."
  ndx = list%find (1)
  write (u, "(A)") "* Test shrinking..."
  call list%shrink_size ()
  write (u, "(A)") "* Test get procedures..."
  call test_get_procedures (list)
  write (u, "(A)") "* Test clearing list..."
  call list%clear ()
  write (u, "(A)") "* Test (more complicated) combinations:"
  write (u, "(A)") "* Test growing (necessary) during adding..."
  call test_grow_and_add (list)
  write (u, "(A)") "* Test adding random data and sorting..."
  call test_sort (list)
  write (u, "(A)") "* Test finding (before sorted)..."
  call test_find (list)
contains
  subroutine test_get_procedures (list)
    type(array_list_t), intent(in) :: list
    integer :: n
    logical :: flag
    n = list%get(1)
    n = list%get_size ()
    n = list%get_count ()
    flag = list%is_element (1)
  end subroutine test_get_procedures

  subroutine test_grow_and_add (list)
    type(array_list_t), intent(inout) :: list
    integer :: i
    do i = 1, 2 * list%get_size ()
       call list%add (i)
    end do
  end subroutine test_grow_and_add

  subroutine test_get (list)
    class(array_list_t), intent(inout) :: list
    integer :: i, data
    do i = list%get_count (), 1, -1
       data = list%get (i)
       if (data == 0) then
          write (u, "(A,1X,I3)") "INDEX EMPTY", i
       end if
    end do
  end subroutine test_get

  subroutine test_sort (list)
    class(array_list_t), intent(inout) :: list
    call list%add (6)
    call list%add (2)
    call list%add (9)
    call list%add (4)
    call list%add (8)
    call list%add (7)
    call list%sort ()
  end subroutine test_sort

  subroutine test_find (list)
    class(array_list_t), intent(inout) :: list
    write (u, "(A,1X,I3)") " 6 INDEX", list%find (6)
    write (u, "(A,1X,I3)") "-1 INDEX", list%find (-1)
    write (u, "(A,1X,I3)") " 3 INDEX", list%find (3)
    write (u, "(A,1X,I3)") "26 INDEX", list%find (26)
    call list%write (u)
  end subroutine test_find
end subroutine array_list_1


end module array_list_uti
