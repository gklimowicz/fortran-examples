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

submodule (user_files) user_files_s

  use io_units
  use diagnostics
  use analysis

  implicit none

contains

  subroutine file_init (file, name, action, status, position)
    type(file_t), intent(out) :: file
    type(string_t), intent(in) :: name
    character(len=*), intent(in) :: action, status, position
    file%unit = free_unit ()
    file%name = name
    open (unit = file%unit, file = char (file%name), &
          action = action, status = status, position = position)
    select case (action)
    case ("read")
       file%reading = .true.
    case ("write")
       file%writing = .true.
    case ("readwrite")
       file%reading = .true.
       file%writing = .true.
    end select
  end subroutine file_init

  subroutine file_final (file)
    type(file_t), intent(inout) :: file
    close (unit = file%unit)
    file%unit = -1
  end subroutine file_final

  function file_is_open (file, action) result (flag)
    logical :: flag
    type(file_t), intent(in) :: file
    character(*), intent(in) :: action
    select case (action)
    case ("read")
       flag = file%reading
    case ("write")
       flag = file%writing
    case ("readwrite")
       flag = file%reading .and. file%writing
    case default
       call msg_bug ("Checking file '" // char (file%name) &
            // "': illegal action specifier")
    end select
  end function file_is_open

  function file_get_unit (file) result (unit)
    integer :: unit
    type(file_t), intent(in) :: file
    unit = file%unit
  end function file_get_unit

  subroutine file_write_string (file, string, advancing)
    type(file_t), intent(in) :: file
    type(string_t), intent(in), optional :: string
    logical, intent(in), optional :: advancing
    if (file%writing) then
       if (present (string)) then
          if (present (advancing)) then
             if (advancing) then
                write (file%unit, "(A)")  char (string)
             else
                write (file%unit, "(A)", advance="no")  char (string)
             end if
          else
             write (file%unit, "(A)")  char (string)
          end if
       else
          write (file%unit, *)
       end if
    else
       call msg_error ("Writing to file: File '" // char (file%name) &
            // "' is not open for writing.")
    end if
  end subroutine file_write_string

  subroutine file_write_ifile (file, ifile)
    type(file_t), intent(in) :: file
    type(ifile_t), intent(in) :: ifile
    type(line_p) :: line
    call line_init (line, ifile)
    do while (line_is_associated (line))
       call file_write_string (file, line_get_string_advance (line))
    end do
  end subroutine file_write_ifile

  subroutine file_write_analysis (file, tag)
    type(file_t), intent(in) :: file
    type(string_t), intent(in), optional :: tag
    if (file%writing) then
       if (present (tag)) then
          call analysis_write (tag, unit = file%unit)
       else
          call analysis_write (unit = file%unit)
       end if
    else
       call msg_error ("Writing analysis to file: File '" // char (file%name) &
            // "' is not open for writing.")
    end if
  end subroutine file_write_analysis

  module subroutine file_list_final (file_list)
    type(file_list_t), intent(inout) :: file_list
    type(file_t), pointer :: current
    do while (associated (file_list%first))
       current => file_list%first
       file_list%first => current%next
       call file_final (current)
       deallocate (current)
    end do
    file_list%last => null ()
  end subroutine file_list_final

  function file_list_get_file_ptr (file_list, name) result (current)
    type(file_t), pointer :: current
    type(file_list_t), intent(in) :: file_list
    type(string_t), intent(in) :: name
    current => file_list%first
    do while (associated (current))
       if (current%name == name)  return
       current => current%next
    end do
  end function file_list_get_file_ptr

  module function file_list_is_open (file_list, name, action) result (flag)
    logical :: flag
    type(file_list_t), intent(in) :: file_list
    type(string_t), intent(in) :: name
    character(len=*), intent(in) :: action
    type(file_t), pointer :: current
    current => file_list_get_file_ptr (file_list, name)
    if (associated (current)) then
       flag = file_is_open (current, action)
    else
       flag = .false.
    end if
  end function file_list_is_open

  module function file_list_get_unit (file_list, name) result (unit)
    integer :: unit
    type(file_list_t), intent(in) :: file_list
    type(string_t), intent(in) :: name
    type(file_t), pointer :: current
    current => file_list_get_file_ptr (file_list, name)
    if (associated (current)) then
       unit = file_get_unit (current)
    else
       unit = -1
    end if
  end function file_list_get_unit

  module subroutine file_list_open (file_list, name, action, status, position)
    type(file_list_t), intent(inout) :: file_list
    type(string_t), intent(in) :: name
    character(len=*), intent(in) :: action, status, position
    type(file_t), pointer :: current
    if (.not. associated (file_list_get_file_ptr (file_list, name))) then
       allocate (current)
       call msg_message ("Opening file '" // char (name) // "' for output")
       call file_init (current, name, action, status, position)
       if (associated (file_list%last)) then
          file_list%last%next => current
          current%prev => file_list%last
       else
          file_list%first => current
       end if
       file_list%last => current
    else
       call msg_error ("Opening file: File '" // char (name) &
            // "' is already open.")
    end if
  end subroutine file_list_open

  module subroutine file_list_close (file_list, name)
    type(file_list_t), intent(inout) :: file_list
    type(string_t), intent(in) :: name
    type(file_t), pointer :: current
    current => file_list_get_file_ptr (file_list, name)
    if (associated (current)) then
       if (associated (current%prev)) then
          current%prev%next => current%next
       else
          file_list%first => current%next
       end if
       if (associated (current%next)) then
          current%next%prev => current%prev
       else
          file_list%last => current%prev
       end if
       call msg_message ("Closing file '" // char (name) // "' for output")
       call file_final (current)
       deallocate (current)
    else
       call msg_error ("Closing file: File '" // char (name) &
            // "' is not open.")
    end if
  end subroutine file_list_close

  module subroutine file_list_write_string (file_list, name, string, advancing)
    type(file_list_t), intent(in) :: file_list
    type(string_t), intent(in) :: name
    type(string_t), intent(in), optional :: string
    logical, intent(in), optional :: advancing
    type(file_t), pointer :: current
    current => file_list_get_file_ptr (file_list, name)
    if (associated (current)) then
       call file_write_string (current, string, advancing)
    else
       call msg_error ("Writing to file: File '" // char (name) &
            // "'is not open.")
    end if
  end subroutine file_list_write_string

  module subroutine file_list_write_ifile (file_list, name, ifile)
    type(file_list_t), intent(in) :: file_list
    type(string_t), intent(in) :: name
    type(ifile_t), intent(in) :: ifile
    type(file_t), pointer :: current
    current => file_list_get_file_ptr (file_list, name)
    if (associated (current)) then
       call file_write_ifile (current, ifile)
    else
       call msg_error ("Writing to file: File '" // char (name) &
            // "'is not open.")
    end if
  end subroutine file_list_write_ifile

  module subroutine file_list_write_analysis (file_list, name, tag)
    type(file_list_t), intent(in) :: file_list
    type(string_t), intent(in) :: name
    type(string_t), intent(in), optional :: tag
    type(file_t), pointer :: current
    if (name == "") then
       if (present (tag)) then
          call analysis_write (tag)
       else
          call analysis_write
       end if
    else
       current => file_list_get_file_ptr (file_list, name)
       if (associated (current)) then
          call file_write_analysis (current, tag)
       else
          call msg_error ("Writing analysis to file: File '" // char (name) &
               // "' is not open.")
       end if
    end if
  end subroutine file_list_write_analysis


end submodule user_files_s

