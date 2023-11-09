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

module file_registries

  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: file_registry_t

  type :: file_handle_t
     type(string_t) :: file
     integer :: unit = 0
     integer :: refcount = 0
   contains
     procedure :: write => file_handle_write
     procedure :: init => file_handle_init
     procedure :: open => file_handle_open
     procedure :: close => file_handle_close
     procedure :: is_open => file_handle_is_open
     procedure :: get_file => file_handle_get_file
     procedure :: get_unit => file_handle_get_unit
  end type file_handle_t

  type, extends (file_handle_t) :: file_entry_t
     type(file_entry_t), pointer :: prev => null ()
     type(file_entry_t), pointer :: next => null ()
  end type file_entry_t

  type :: file_registry_t
     type(file_entry_t), pointer :: first => null ()
   contains
     procedure :: write => file_registry_write
     procedure :: open => file_registry_open
     procedure :: close => file_registry_close
  end type file_registry_t


  interface
    module subroutine file_handle_write (handle, u, show_unit)
      class(file_handle_t), intent(in) :: handle
      integer, intent(in) :: u
      logical, intent(in), optional :: show_unit
    end subroutine file_handle_write
    module subroutine file_handle_init (handle, file)
      class(file_handle_t), intent(out) :: handle
      type(string_t), intent(in) :: file
    end subroutine file_handle_init
    module subroutine file_handle_open (handle)
      class(file_handle_t), intent(inout) :: handle
    end subroutine file_handle_open
    module subroutine file_handle_close (handle)
      class(file_handle_t), intent(inout) :: handle
    end subroutine file_handle_close
    module function file_handle_is_open (handle) result (flag)
      class(file_handle_t), intent(in) :: handle
      logical :: flag
    end function file_handle_is_open
    module function file_handle_get_file (handle) result (file)
      class(file_handle_t), intent(in) :: handle
      type(string_t) :: file
    end function file_handle_get_file
    module function file_handle_get_unit (handle) result (unit)
      class(file_handle_t), intent(in) :: handle
      integer :: unit
    end function file_handle_get_unit
    module subroutine file_registry_write (registry, unit, show_unit)
      class(file_registry_t), intent(in) :: registry
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_unit
    end subroutine file_registry_write
    module subroutine file_registry_open (registry, file, unit)
      class(file_registry_t), intent(inout) :: registry
      type(string_t), intent(in) :: file
      integer, intent(out), optional :: unit
    end subroutine file_registry_open
  module subroutine file_registry_close (registry, file)
    class(file_registry_t), intent(inout) :: registry
    type(string_t), intent(in) :: file
  end subroutine file_registry_close
  end interface

end module file_registries
