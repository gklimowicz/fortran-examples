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

module prclib_stacks_uti

  use iso_varying_string, string_t => varying_string

  use prclib_stacks

  implicit none
  private

  public :: prclib_stacks_1
  public :: prclib_stacks_2

contains

  subroutine prclib_stacks_1 (u)
    integer, intent(in) :: u
    type(prclib_stack_t) :: stack

    write (u, "(A)")  "* Test output: prclib_stacks_1"
    write (u, "(A)")  "*   Purpose: display an empty process library stack"
    write (u, "(A)")

    call stack%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prclib_stacks_1"

  end subroutine prclib_stacks_1

  subroutine prclib_stacks_2 (u)
    integer, intent(in) :: u
    type(prclib_stack_t) :: stack
    type(prclib_entry_t), pointer :: lib

    write (u, "(A)")  "* Test output: prclib_stacks_2"
    write (u, "(A)")  "*   Purpose: fill a process library stack"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize two (empty) libraries &
         &and push them on the stack"
    write (u, "(A)")

    allocate (lib)
    call lib%init (var_str ("lib1"))
    call stack%push (lib)

    allocate (lib)
    call lib%init (var_str ("lib2"))
    call stack%push (lib)

    call stack%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call stack%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: prclib_stacks_2"

  end subroutine prclib_stacks_2


end module prclib_stacks_uti
