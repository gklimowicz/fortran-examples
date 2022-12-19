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

module array_list
  use kinds, only: default

  implicit none
  private

  public :: array_list_t

  integer, parameter :: ARRAY_LIST_START_SIZE = 10
  real(default), parameter :: ARRAY_LIST_GROW_FACTOR = 1.5_default, &
       ARRAY_LIST_SHRINK_THRESHOLD = 0.3_default


  type :: array_list_t
     private
     integer, dimension(:), allocatable :: array
     !! Track the index to *current* item, to be stored.
     !! Must fulfill: 0 <= count <= size.
     integer :: count = 0
     !! size \in N.
     integer :: size = 0
   contains
     procedure :: write => array_list_write
     procedure :: init => array_list_init
     procedure :: get => array_list_get
     procedure :: get_count => array_list_get_count
     procedure :: get_size => array_list_get_size
     procedure :: is_full => array_list_is_full
     procedure :: is_empty => array_list_is_empty
     procedure :: is_index => array_list_is_index
     procedure :: clear => array_list_clear
     procedure :: add => array_list_add
     procedure :: grow_size => array_list_grow_size
     procedure :: shrink_size => array_list_shrink_size
     procedure :: reverse_order => array_list_reverse_order
     procedure :: sort => array_list_sort
     procedure :: is_element => array_list_is_element
     procedure :: find => array_list_find
     procedure :: add_at => array_list_add_at
     procedure :: remove => array_list_remove
     procedure :: remove_at => array_list_remove_at
  end type array_list_t


  interface
    module subroutine array_list_write (list, unit)
      class(array_list_t), intent(in) :: list
      integer, intent(in), optional :: unit
    end subroutine array_list_write
    module subroutine array_list_init (list)
      class(array_list_t), intent(out) :: list
    end subroutine array_list_init
    elemental module function array_list_get (list, index) result (data)
      class(array_list_t), intent(in) :: list
      integer, intent(in) :: index
      integer :: data
    end function array_list_get
    pure module function array_list_get_count (list) result (count)
      class(array_list_t), intent(in) :: list
      integer :: count
    end function array_list_get_count
    pure module function array_list_get_size (list) result (size)
      class(array_list_t), intent(in) :: list
      integer :: size
    end function array_list_get_size
    pure module function array_list_is_full (list) result (flag)
      class(array_list_t), intent(in) :: list
      logical :: flag
    end function array_list_is_full
    pure module function array_list_is_empty (list) result (flag)
      class(array_list_t), intent(in) :: list
      logical :: flag
    end function array_list_is_empty
    pure module function array_list_is_index (list, index) result (flag)
      class(array_list_t), intent(in) :: list
      integer, intent(in) :: index
      logical :: flag
    end function array_list_is_index
    module subroutine array_list_clear (list)
      class(array_list_t), intent(inout) :: list
    end subroutine array_list_clear
    module subroutine array_list_add (list, data)
      class(array_list_t), intent(inout) :: list
      integer, intent(in) :: data
    end subroutine array_list_add
    module subroutine array_list_grow_size (list)
      class(array_list_t), intent(inout) :: list
    end subroutine array_list_grow_size
    module subroutine array_list_shrink_size (list)
      class(array_list_t), intent(inout) :: list
      integer, dimension(:), allocatable :: array
    end subroutine array_list_shrink_size
    module subroutine array_list_reverse_order (list)
      class(array_list_t), intent(inout) :: list
    end subroutine array_list_reverse_order
    pure module subroutine array_list_sort (list)
      class(array_list_t), intent(inout) :: list
    end subroutine array_list_sort
    pure module function array_list_is_element (list, data) result (flag)
      class(array_list_t), intent(in) :: list
      integer, intent(in) :: data
      logical :: flag
    end function array_list_is_element
    module function array_list_find (list, data) result (index)
      class(array_list_t), intent(inout) :: list
      integer, intent(in) :: data
      integer :: index
    end function array_list_find
    module subroutine array_list_add_at (list, index, data)
      class(array_list_t), intent(inout) :: list
      integer, intent(in) :: index
      integer, intent(in) :: data
    end subroutine array_list_add_at
    module function array_list_remove (list) result (data)
      class(array_list_t), intent(inout) :: list
      integer :: data
    end function array_list_remove
    module function array_list_remove_at (list, index) result (data)
      class(array_list_t), intent(inout) :: list
      integer, intent(in) :: index
      integer :: data
    end function array_list_remove_at
  end interface

end module array_list
