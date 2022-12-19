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

module process_stacks

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use variables
  use process

  implicit none
  private

  public :: process_entry_t
  public :: process_stack_t

  type, extends (process_t) :: process_entry_t
     type(process_entry_t), pointer :: next => null ()
  end type process_entry_t

  type :: process_stack_t
     integer :: n = 0
     type(process_entry_t), pointer :: first => null ()
     type(var_list_t), pointer :: var_list => null ()
     type(process_stack_t), pointer :: next => null ()
   contains
     procedure :: clear => process_stack_clear
     procedure :: final => process_stack_final
     procedure :: write => process_stack_write
     procedure :: write_var_list => process_stack_write_var_list
     procedure :: show => process_stack_show
     procedure :: link => process_stack_link
     procedure :: init_var_list => process_stack_init_var_list
     procedure :: link_var_list => process_stack_link_var_list
     procedure :: push => process_stack_push
     procedure :: pop_last => process_stack_pop_last
     procedure :: init_result_vars => process_stack_init_result_vars
     procedure :: fill_result_vars => process_stack_fill_result_vars
     procedure :: update_result_vars => process_stack_update_result_vars
     procedure :: exists => process_stack_exists
     procedure :: get_process_ptr => process_stack_get_process_ptr
  end type process_stack_t


  interface
    module subroutine process_stack_clear (stack)
      class(process_stack_t), intent(inout) :: stack
    end subroutine process_stack_clear
    module subroutine process_stack_final (object)
      class(process_stack_t), intent(inout) :: object
    end subroutine process_stack_final
    recursive module subroutine process_stack_write (object, unit, pacify)
      class(process_stack_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: pacify
    end subroutine process_stack_write
    module subroutine process_stack_write_var_list (object, unit)
      class(process_stack_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine process_stack_write_var_list
    recursive module subroutine process_stack_show (object, unit, fifo)
      class(process_stack_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: fifo
    end subroutine process_stack_show
    module subroutine process_stack_link (local_stack, global_stack)
      class(process_stack_t), intent(inout) :: local_stack
      type(process_stack_t), intent(in), target :: global_stack
    end subroutine process_stack_link
    module subroutine process_stack_init_var_list (stack, var_list)
      class(process_stack_t), intent(inout) :: stack
      type(var_list_t), intent(inout), optional :: var_list
    end subroutine process_stack_init_var_list
    module subroutine process_stack_link_var_list (stack, var_list)
      class(process_stack_t), intent(inout) :: stack
      type(var_list_t), intent(in), target :: var_list
    end subroutine process_stack_link_var_list
    module subroutine process_stack_push (stack, process)
      class(process_stack_t), intent(inout) :: stack
      type(process_entry_t), intent(inout), pointer :: process
    end subroutine process_stack_push
    module subroutine process_stack_pop_last (stack, process)
      class(process_stack_t), intent(inout) :: stack
      type(process_entry_t), intent(inout), pointer :: process
    end subroutine process_stack_pop_last
    module subroutine process_stack_init_result_vars (stack, id)
      class(process_stack_t), intent(inout) :: stack
      type(string_t), intent(in) :: id
    end subroutine process_stack_init_result_vars
    module subroutine process_stack_fill_result_vars (stack, id)
      class(process_stack_t), intent(inout) :: stack
      type(string_t), intent(in) :: id
    end subroutine process_stack_fill_result_vars
    module subroutine process_stack_update_result_vars &
         (stack, id, var_list_local)
      class(process_stack_t), intent(inout) :: stack
      type(string_t), intent(in) :: id
      type(var_list_t), intent(inout) :: var_list_local
    end subroutine process_stack_update_result_vars
    module function process_stack_exists (stack, id) result (flag)
      class(process_stack_t), intent(in) :: stack
      type(string_t), intent(in) :: id
      logical :: flag
    end function process_stack_exists
    recursive module function process_stack_get_process_ptr &
         (stack, id) result (ptr)
      class(process_stack_t), intent(in) :: stack
      type(string_t), intent(in) :: id
      type(process_t), pointer :: ptr
    end function process_stack_get_process_ptr
  end interface

end module process_stacks
