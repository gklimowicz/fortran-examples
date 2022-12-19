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

submodule (ifiles) ifiles_s

  use system_defs, only: EOF
  use io_units

  implicit none

contains

  subroutine line_entry_create (line, string)
    type(line_entry_t), pointer :: line
    type(string_t), intent(in) :: string
    allocate (line)
    line%string = string
  end subroutine line_entry_create

  subroutine line_entry_destroy (line)
    type(line_entry_t), pointer :: line
    deallocate (line)
  end subroutine line_entry_destroy

  module subroutine ifile_clear (ifile)
    class(ifile_t), intent(inout) :: ifile
    type(line_entry_t), pointer :: current
    do while (associated (ifile%first))
       current => ifile%first
       ifile%first => current%next
       call line_entry_destroy (current)
    end do
    nullify (ifile%last)
    ifile%n_lines = 0
  end subroutine ifile_clear

  module subroutine ifile_read_from_string (ifile, string)
    type(ifile_t), intent(inout) :: ifile
    type(string_t), intent(in) :: string
    call ifile_clear (ifile)
    call ifile_append (ifile, string)
  end subroutine ifile_read_from_string

  module subroutine ifile_read_from_char (ifile, char)
    type(ifile_t), intent(inout) :: ifile
    character(*), intent(in) :: char
    call ifile_clear (ifile)
    call ifile_append (ifile, char)
  end subroutine ifile_read_from_char

  module subroutine ifile_read_from_char_array (ifile, char)
    type(ifile_t), intent(inout) :: ifile
    character(*), dimension(:), intent(in) :: char
    call ifile_clear (ifile)
    call ifile_append (ifile, char)
  end subroutine ifile_read_from_char_array

  module subroutine ifile_read_from_unit (ifile, unit, iostat)
    type(ifile_t), intent(inout) :: ifile
    integer, intent(in) :: unit
    integer, intent(out), optional :: iostat
    call ifile_clear (ifile)
    call ifile_append (ifile, unit, iostat)
  end subroutine ifile_read_from_unit

  module subroutine ifile_read_from_ifile (ifile, ifile_in)
    type(ifile_t), intent(inout) :: ifile
    type(ifile_t), intent(in) :: ifile_in
    call ifile_clear (ifile)
    call ifile_append (ifile, ifile_in)
  end subroutine ifile_read_from_ifile

  module subroutine ifile_append_from_string (ifile, string)
    class(ifile_t), intent(inout) :: ifile
    type(string_t), intent(in) :: string
    type(line_entry_t), pointer :: current
    call line_entry_create (current, string)
    current%index = ifile%n_lines + 1
    if (associated (ifile%last)) then
       current%previous => ifile%last
       ifile%last%next => current
    else
       ifile%first => current
    end if
    ifile%last => current
    ifile%n_lines = current%index
  end subroutine ifile_append_from_string

  module subroutine ifile_append_from_char (ifile, char)
    class(ifile_t), intent(inout) :: ifile
    character(*), intent(in) :: char
    call ifile_append_from_string (ifile, var_str (trim (char)))
  end subroutine ifile_append_from_char

  module subroutine ifile_append_from_char_array (ifile, char)
    class(ifile_t), intent(inout) :: ifile
    character(*), dimension(:), intent(in) :: char
    integer :: i
    do i = 1, size (char)
       call ifile_append_from_string (ifile, var_str (trim (char(i))))
    end do
  end subroutine ifile_append_from_char_array

  module subroutine ifile_append_from_unit (ifile, unit, iostat)
    class(ifile_t), intent(inout) :: ifile
    integer, intent(in) :: unit
    integer, intent(out), optional :: iostat
    type(string_t) :: buffer
    integer :: ios
    ios = 0
    READ_LOOP: do
       call get (unit, buffer, iostat = ios)
       if (ios == EOF .or. ios > 0)  exit READ_LOOP
       call ifile_append_from_string (ifile, buffer)
    end do READ_LOOP
    if (present (iostat)) then
       iostat = ios
    else if (ios > 0) then
       call get (unit, buffer)  ! trigger error again
    end if
  end subroutine ifile_append_from_unit

  module subroutine ifile_append_from_ifile (ifile, ifile_in)
    class(ifile_t), intent(inout) :: ifile
    type(ifile_t), intent(in) :: ifile_in
    type(line_entry_t), pointer :: current
    current => ifile_in%first
    do while (associated (current))
       call ifile_append_from_string (ifile, current%string)
       current => current%next
    end do
  end subroutine ifile_append_from_ifile

  module subroutine ifile_write (ifile, unit, iostat)
    type(ifile_t), intent(in) :: ifile
    integer, intent(in), optional :: unit
    integer, intent(out), optional :: iostat
    integer :: u
    type(line_entry_t), pointer :: current
    u = given_output_unit (unit);  if (u < 0)  return
    current => ifile%first
    do while (associated (current))
       call put_line (u, current%string, iostat)
       current => current%next
    end do
  end subroutine ifile_write

  module subroutine ifile_to_string_array (ifile, string)
    type(ifile_t), intent(in) :: ifile
    type(string_t), dimension(:), intent(inout), allocatable :: string
    type(line_entry_t), pointer :: current
    integer :: i
    allocate (string (ifile_get_length (ifile)))
    current => ifile%first
    do i = 1, ifile_get_length (ifile)
       string(i) = current%string
       current => current%next
    end do
  end subroutine ifile_to_string_array

  module function ifile_get_length (ifile) result (length)
    integer :: length
    type(ifile_t), intent(in) :: ifile
    length = ifile%n_lines
  end function ifile_get_length

  module subroutine line_init (line, ifile, back)
    type(line_p), intent(inout) :: line
    type(ifile_t), intent(in) :: ifile
    logical, intent(in), optional :: back
    if (present (back)) then
       if (back) then
          line%p => ifile%last
       else
          line%p => ifile%first
       end if
    else
       line%p => ifile%first
    end if
  end subroutine line_init

  module subroutine line_final (line)
    type(line_p), intent(inout) :: line
    nullify (line%p)
  end subroutine line_final

  module subroutine line_advance (line)
    type(line_p), intent(inout) :: line
    if (associated (line%p))  line%p => line%p%next
  end subroutine line_advance

  module subroutine line_backspace (line)
    type(line_p), intent(inout) :: line
    if (associated (line%p))  line%p => line%p%previous
  end subroutine line_backspace

  module function line_is_associated (line) result (ok)
    logical :: ok
    type(line_p), intent(in) :: line
    ok = associated (line%p)
  end function line_is_associated

  module function line_get_string (line) result (string)
    type(string_t) :: string
    type(line_p), intent(in) :: line
    if (associated (line%p)) then
       string = line%p%string
    else
       string = ""
    end if
  end function line_get_string

  module function line_get_string_advance (line) result (string)
    type(string_t) :: string
    type(line_p), intent(inout) :: line
    if (associated (line%p)) then
       string = line%p%string
       call line_advance (line)
    else
       string = ""
    end if
  end function line_get_string_advance

  module function line_get_index (line) result (index)
    integer :: index
    type(line_p), intent(in) :: line
    if (associated (line%p)) then
       index = line%p%index
    else
       index = 0
    end if
  end function line_get_index

  module function line_get_length (line) result (length)
    integer :: length
    type(line_p), intent(in) :: line
    if (associated (line%p)) then
       length = len (line%p%string)
    else
       length = 0
    end if
  end function line_get_length


end submodule ifiles_s

