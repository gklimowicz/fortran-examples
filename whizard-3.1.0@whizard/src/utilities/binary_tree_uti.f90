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

module binary_tree_uti

  use binary_tree

  implicit none
  private

  type :: btree_obj_t
     integer :: i = 0
  end type btree_obj_t

public :: binary_tree_1

 contains

  subroutine binary_tree_1 (u)
    integer, intent(in) :: u
    integer, dimension(10) :: ndx = [1, 2, 5, 7, 19, 23, 97, -1, -6, 0]
    class(*), pointer :: obj
    type(binary_tree_t) :: btree
    type(binary_tree_iterator_t) :: iterator
    integer :: i, key
    write (u, "(A)") "* Test output: Binary tree"
    write (u, "(A)") "*   Purpose: test interface and implementation of binary tree " // &
         "and its iterator using polymorph objects."
    write (u, "(A)")

    write (u, "(A)") "* Insert fixed number of objects into tree..."
    do i = 1, size (ndx)
       call allocate_obj (i, obj)
       call btree%insert (ndx(i), obj)
    end do

    write (u, "(A)") "* Search for all added objects in tree..."
    do i = size (ndx), 1, -1
       write (u, "(A,1X,I3,1X,L1)") "- Has key", ndx(i), btree%has_key (ndx(i))
       call btree%search (ndx(i), obj)
       select type (obj)
       type is (btree_obj_t)
          write (u, "(2(A,1X,I3,1X))") "- NDX", ndx(i), "OBJ", obj%i
       end select
    end do

    write (u, "(A)") "* Output binary tree in preorder..."
    call btree%write (u)

    write (u, "(A)") "* Clear binary tree..."
    call btree%clear ()
    call btree%write (u)

    write (u, "(A)") "* Insert fixed number of object into tree (reversed order)..."
    do i = size (ndx), 1, -1
       call allocate_obj (i, obj)
       call btree%insert (ndx(i), obj)
    end do

    write (u, "(A)") "* Iterate over binary tree..."
    call iterator%init (btree)
    do while (iterator%is_iterable ())
       call iterator%next (key)
       call btree%search (key, obj)
       select type (obj)
       type is (btree_obj_t)
          write (u, "(2(A,1X,I3,1X))") "- KEY", key, "OBJ", obj%i
       end select
    end do

    write (u, "(A)") "* Search for a non-existing key..."
    write (u, "(A,1X,I3,1X,L1)") "- Has key", 123, btree%has_key (123)
    call btree%search (123, obj)
    write (u, "(A,1X,L1)") "- Object found", associated (obj)

    !! Do not test against a duplicate entry as the it will forcibly stop the program.
  contains
    subroutine allocate_obj (num, obj)
      integer, intent(in) :: num
      class(*), pointer, intent(out) :: obj
      allocate (btree_obj_t :: obj)
      select type (obj)
      type is (btree_obj_t)
         obj%i = num
      end select
    end subroutine allocate_obj
  end subroutine binary_tree_1


end module binary_tree_uti
