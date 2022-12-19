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
module binary_tree

  implicit none
  private

  public :: binary_tree_iterator_t 
  public :: binary_tree_t
  
  type :: binary_tree_iterator_t
     integer, dimension(:), allocatable :: key
     integer :: current
     !! current \in {1, N}.
   contains
     procedure :: init => binary_tree_iterator_init
     procedure :: is_iterable => binary_tree_iterator_is_iterable
     procedure :: next => binary_tree_iterator_next
  end type binary_tree_iterator_t

  type :: binary_tree_node_t
     integer :: height = 0
     type(binary_tree_node_t), pointer :: left => null ()
     type(binary_tree_node_t), pointer :: right => null ()
     !!
     integer :: key = 0
     class(*), pointer :: obj => null ()
   contains
     procedure :: init => binary_tree_node_init
     procedure :: write => binary_tree_node_write
     procedure :: get_balance => binary_tree_node_get_balance
     procedure :: increment_height => binary_tree_node_increment_height
     final :: binary_tree_node_final
  end type binary_tree_node_t

  type :: binary_tree_t
     integer :: n_elements = 0
     type(binary_tree_node_t), pointer :: root => null ()
   contains
     procedure :: write => binary_tree_write
     final :: binary_tree_final
     procedure :: clear => binary_tree_clear
     procedure :: get_n_elements => binary_tree_get_n_elements
     procedure :: insert => binary_tree_insert
     procedure, private :: insert_node => binary_tree_insert_node
     procedure, private :: balance => binary_tree_balance
     procedure :: search => binary_tree_search
     procedure :: has_key => binary_tree_has_key
     procedure, private :: rotate_right => binary_tree_rotate_right
     procedure, private :: rotate_left => binary_tree_rotate_left
  end type binary_tree_t


  interface
    module subroutine binary_tree_iterator_init (iterator, btree)
      class(binary_tree_iterator_t), intent(inout) :: iterator
      type(binary_tree_t), target :: btree
    end subroutine binary_tree_iterator_init
    module function binary_tree_iterator_is_iterable (iterator) result (flag)
      class(binary_tree_iterator_t), intent(in) :: iterator
      logical :: flag
    end function binary_tree_iterator_is_iterable
    module subroutine binary_tree_iterator_next (iterator, key)
      class(binary_tree_iterator_t), intent(inout) :: iterator
      integer, intent(out) :: key
    end subroutine binary_tree_iterator_next
    module subroutine binary_tree_node_init (btree_node, key, obj)
      class(binary_tree_node_t), intent(inout) :: btree_node
      integer, intent(in) :: key
      class(*), pointer :: obj
    end subroutine binary_tree_node_init
    recursive module subroutine binary_tree_node_write &
         (btree_node, unit, level, mode)
      class(binary_tree_node_t), intent(in) :: btree_node
      integer, intent(in) :: unit
      integer, intent(in) :: level
      character(len=*), intent(in) :: mode
    end subroutine binary_tree_node_write
    module function binary_tree_node_get_balance (btree_node) result (balance)
      class(binary_tree_node_t), intent(in) :: btree_node
      integer :: balance
    end function binary_tree_node_get_balance
    module subroutine binary_tree_node_increment_height (btree_node)
      class(binary_tree_node_t), intent(inout) :: btree_node
    end subroutine binary_tree_node_increment_height
    !!! !!! NAG 7 compiler bug with finalizers and unlimited polymorphism
    !!! module subroutine binary_tree_node_final (btree_node)
    !!!   type(binary_tree_node_t), intent(inout) :: btree_node
    !!! end subroutine binary_tree_node_final
    module subroutine binary_tree_write (btree, unit)
      class(binary_tree_t), intent(in) :: btree
      integer, intent(in), optional :: unit
    end subroutine binary_tree_write
    !!! !!! NAG 7 compiler bug with finalizers and unlimited polymorphism
    !!! module subroutine binary_tree_final (btree)
    !!!   type(binary_tree_t), intent(inout) :: btree
    !!! end subroutine binary_tree_final
    module subroutine binary_tree_clear (btree)
      class(binary_tree_t), intent(inout) :: btree
    end subroutine binary_tree_clear
    module function binary_tree_get_n_elements (btree) result (n)
      class(binary_tree_t), intent(in) :: btree
      integer :: n
    end function binary_tree_get_n_elements
    module subroutine binary_tree_insert (btree, key, obj)
      class(binary_tree_t), intent(inout) :: btree
      integer, intent(in) :: key
      class(*), pointer, intent(in) :: obj
    end subroutine binary_tree_insert
    recursive module subroutine binary_tree_insert_node (btree, parent, node)
      class(binary_tree_t), intent(in) :: btree
      type(binary_tree_node_t), intent(inout), pointer :: parent
      type(binary_tree_node_t), intent(in), pointer :: node
    end subroutine binary_tree_insert_node
    module subroutine binary_tree_balance (btree, subtree, key)
      class(binary_tree_t), intent(in) :: btree
      type(binary_tree_node_t), intent(inout), pointer :: subtree
      integer, intent(in) :: key
    end subroutine binary_tree_balance
    module subroutine binary_tree_search (btree, key, obj)
      class(binary_tree_t), intent(in) :: btree
      integer, intent(in) :: key
      class(*), pointer, intent(out) :: obj
    end subroutine binary_tree_search
    module function binary_tree_has_key (btree, key) result (flag)
      class(binary_tree_t), intent(in) :: btree
      integer, intent(in) :: key
      logical :: flag
    end function binary_tree_has_key
    module subroutine binary_tree_rotate_right (btree, root, new_root)
      class(binary_tree_t), intent(in) :: btree
      type(binary_tree_node_t), pointer, intent(inout) :: root
      type(binary_tree_node_t), pointer, intent(out) :: new_root
    end subroutine binary_tree_rotate_right
    module subroutine binary_tree_rotate_left (btree, root, new_root)
      class(binary_tree_t), intent(in) :: btree
      type(binary_tree_node_t), pointer, intent(inout) :: root
      type(binary_tree_node_t), pointer, intent(out) :: new_root
    end subroutine binary_tree_rotate_left
  end interface

contains

  recursive subroutine binary_tree_node_final (btree_node)
    type(binary_tree_node_t), intent(inout) :: btree_node
    if (associated (btree_node%left)) deallocate (btree_node%left)
    if (associated (btree_node%right)) deallocate (btree_node%right)
    deallocate (btree_node%obj)
  end subroutine binary_tree_node_final

  subroutine binary_tree_final (btree)
    type(binary_tree_t), intent(inout) :: btree
    btree%n_elements = 0
    if (associated (btree%root)) then
       deallocate (btree%root)
    end if
  end subroutine binary_tree_final

  
end module binary_tree
