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

module prclib_stacks

  use iso_varying_string, string_t => varying_string
  use process_libraries

  implicit none
  private

  public :: prclib_entry_t
  public :: prclib_stack_t

  type, extends (process_library_t) :: prclib_entry_t
     type(prclib_entry_t), pointer :: next => null ()
  end type prclib_entry_t

  type :: prclib_stack_t
     integer :: n = 0
     type(prclib_entry_t), pointer :: first => null ()
   contains
     procedure :: final => prclib_stack_final
     procedure :: write => prclib_stack_write
     procedure :: push => prclib_stack_push
     procedure :: get_first_ptr => prclib_stack_get_first_ptr
     procedure :: get_names => prclib_stack_get_names
     procedure :: get_library_ptr => prclib_stack_get_library_ptr
  end type prclib_stack_t


  interface
    module subroutine prclib_stack_final (object)
      class(prclib_stack_t), intent(inout) :: object
    end subroutine prclib_stack_final
    module subroutine prclib_stack_write (object, unit, libpath)
      class(prclib_stack_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: libpath
    end subroutine prclib_stack_write
    module subroutine prclib_stack_push (stack, lib)
      class(prclib_stack_t), intent(inout) :: stack
      type(prclib_entry_t), intent(inout), pointer :: lib
    end subroutine prclib_stack_push
    module function prclib_stack_get_first_ptr (stack) result (ptr)
      class(prclib_stack_t), intent(in) :: stack
      type(process_library_t), pointer :: ptr
    end function prclib_stack_get_first_ptr
    module subroutine prclib_stack_get_names (stack, libname)
      class(prclib_stack_t), intent(in) :: stack
      type(string_t), dimension(:), allocatable, intent(out) :: libname
    end subroutine prclib_stack_get_names
    module function prclib_stack_get_library_ptr (stack, libname) result (ptr)
      class(prclib_stack_t), intent(in) :: stack
      type(string_t), intent(in) :: libname
      type(process_library_t), pointer :: ptr
    end function prclib_stack_get_library_ptr
  end interface

end module prclib_stacks
