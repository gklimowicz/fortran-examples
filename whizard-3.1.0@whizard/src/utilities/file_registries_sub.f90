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

submodule (file_registries) file_registries_s

  use io_units

  implicit none

contains

  module subroutine file_handle_write (handle, u, show_unit)
    class(file_handle_t), intent(in) :: handle
    integer, intent(in) :: u
    logical, intent(in), optional :: show_unit
    logical :: show_u
    show_u = .false.;  if (present (show_unit))  show_u = show_unit
    if (show_u) then
       write (u, "(3x,A,1x,I0,1x,'(',I0,')')")  &
            char (handle%file), handle%unit, handle%refcount
    else
       write (u, "(3x,A,1x,'(',I0,')')")  &
            char (handle%file), handle%refcount
    end if
  end subroutine file_handle_write

  module subroutine file_handle_init (handle, file)
    class(file_handle_t), intent(out) :: handle
    type(string_t), intent(in) :: file
    handle%file = file
  end subroutine file_handle_init

  module subroutine file_handle_open (handle)
    class(file_handle_t), intent(inout) :: handle
    if (handle%refcount == 0) then
       handle%unit = free_unit ()
       open (unit = handle%unit, file = char (handle%file), action = "read", &
            status = "old")
    end if
    handle%refcount = handle%refcount + 1
  end subroutine file_handle_open

  module subroutine file_handle_close (handle)
    class(file_handle_t), intent(inout) :: handle
    handle%refcount = handle%refcount - 1
    if (handle%refcount == 0) then
       close (handle%unit)
       handle%unit = 0
    end if
  end subroutine file_handle_close

  module function file_handle_is_open (handle) result (flag)
    class(file_handle_t), intent(in) :: handle
    logical :: flag
    flag = handle%unit /= 0
  end function file_handle_is_open

  module function file_handle_get_file (handle) result (file)
    class(file_handle_t), intent(in) :: handle
    type(string_t) :: file
    file = handle%file
  end function file_handle_get_file

  module function file_handle_get_unit (handle) result (unit)
    class(file_handle_t), intent(in) :: handle
    integer :: unit
    unit = handle%unit
  end function file_handle_get_unit

  module subroutine file_registry_write (registry, unit, show_unit)
    class(file_registry_t), intent(in) :: registry
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_unit
    type(file_entry_t), pointer :: entry
    integer :: u
    u = given_output_unit (unit)
    if (associated (registry%first)) then
       write (u, "(1x,A)")  "File registry:"
       entry => registry%first
       do while (associated (entry))
          call entry%write (u, show_unit)
          entry => entry%next
       end do
    else
       write (u, "(1x,A)")  "File registry: [empty]"
    end if
  end subroutine file_registry_write

  module subroutine file_registry_open (registry, file, unit)
    class(file_registry_t), intent(inout) :: registry
    type(string_t), intent(in) :: file
    integer, intent(out), optional :: unit
    type(file_entry_t), pointer :: entry
    entry => registry%first
    FIND_ENTRY: do while (associated (entry))
       if (entry%get_file () == file)  exit FIND_ENTRY
       entry => entry%next
    end do FIND_ENTRY
    if (.not. associated (entry)) then
       allocate (entry)
       call entry%init (file)
       if (associated (registry%first)) then
          registry%first%prev => entry
          entry%next => registry%first
       end if
       registry%first => entry
    end if
    call entry%open ()
    if (present (unit))  unit = entry%get_unit ()
  end subroutine file_registry_open

  module subroutine file_registry_close (registry, file)
    class(file_registry_t), intent(inout) :: registry
    type(string_t), intent(in) :: file
    type(file_entry_t), pointer :: entry
    entry => registry%first
    FIND_ENTRY: do while (associated (entry))
       if (entry%get_file () == file)  exit FIND_ENTRY
       entry => entry%next
    end do FIND_ENTRY
    if (associated (entry)) then
       call entry%close ()
       if (.not. entry%is_open ()) then
          if (associated (entry%prev)) then
             entry%prev%next => entry%next
          else
             registry%first => entry%next
          end if
          if (associated (entry%next)) then
             entry%next%prev => entry%prev
          end if
          deallocate (entry)
       end if
    end if
  end subroutine file_registry_close


end submodule file_registries_s
