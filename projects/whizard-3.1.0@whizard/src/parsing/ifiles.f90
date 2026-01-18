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

module ifiles

  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: ifile_t
  public :: ifile_clear
  public :: ifile_final
  public :: ifile_read
  public :: ifile_append
  public :: ifile_write
  public :: ifile_to_string_array
  public :: ifile_get_length
  public :: line_p
  public :: line_init
  public :: line_final
  public :: line_advance
  public :: line_backspace
  public :: line_is_associated
  public :: line_get_string
  public :: line_get_string_advance
  public :: line_get_index
  public :: line_get_length

  type :: line_entry_t
     private
     type(line_entry_t), pointer :: previous => null ()
     type(line_entry_t), pointer :: next => null ()
     type(string_t) :: string
     integer :: index
  end type line_entry_t

  type :: ifile_t
     private
     type(line_entry_t), pointer :: first => null ()
     type(line_entry_t), pointer :: last => null ()
     integer :: n_lines = 0
   contains
     procedure :: final => ifile_clear
     generic :: append => &
          ifile_append_from_char
     procedure, private :: ifile_append_from_char
  end type ifile_t

  type :: line_p
     private
     type(line_entry_t), pointer :: p => null ()
  end type line_p


  interface ifile_final
     module procedure ifile_clear
  end interface
  interface ifile_read
     module procedure ifile_read_from_string
     module procedure ifile_read_from_char
     module procedure ifile_read_from_unit
     module procedure ifile_read_from_char_array
     module procedure ifile_read_from_ifile
  end interface
  interface ifile_append
     module procedure ifile_append_from_string
     module procedure ifile_append_from_char
     module procedure ifile_append_from_unit
     module procedure ifile_append_from_char_array
     module procedure ifile_append_from_ifile
  end interface

  interface
    module subroutine ifile_clear (ifile)
      class(ifile_t), intent(inout) :: ifile
    end subroutine ifile_clear
    module subroutine ifile_read_from_string (ifile, string)
      type(ifile_t), intent(inout) :: ifile
      type(string_t), intent(in) :: string
    end subroutine ifile_read_from_string
    module subroutine ifile_read_from_char (ifile, char)
      type(ifile_t), intent(inout) :: ifile
      character(*), intent(in) :: char
    end subroutine ifile_read_from_char
    module subroutine ifile_read_from_char_array (ifile, char)
      type(ifile_t), intent(inout) :: ifile
      character(*), dimension(:), intent(in) :: char
    end subroutine ifile_read_from_char_array
    module subroutine ifile_read_from_unit (ifile, unit, iostat)
      type(ifile_t), intent(inout) :: ifile
      integer, intent(in) :: unit
      integer, intent(out), optional :: iostat
    end subroutine ifile_read_from_unit
    module subroutine ifile_read_from_ifile (ifile, ifile_in)
      type(ifile_t), intent(inout) :: ifile
      type(ifile_t), intent(in) :: ifile_in
    end subroutine ifile_read_from_ifile
    module subroutine ifile_append_from_string (ifile, string)
      class(ifile_t), intent(inout) :: ifile
      type(string_t), intent(in) :: string
    end subroutine ifile_append_from_string
    module subroutine ifile_append_from_char (ifile, char)
      class(ifile_t), intent(inout) :: ifile
      character(*), intent(in) :: char
    end subroutine ifile_append_from_char
    module subroutine ifile_append_from_char_array (ifile, char)
      class(ifile_t), intent(inout) :: ifile
      character(*), dimension(:), intent(in) :: char
    end subroutine ifile_append_from_char_array
    module subroutine ifile_append_from_ifile (ifile, ifile_in)
      class(ifile_t), intent(inout) :: ifile
      type(ifile_t), intent(in) :: ifile_in
    end subroutine ifile_append_from_ifile
    module subroutine ifile_append_from_unit (ifile, unit, iostat)
      class(ifile_t), intent(inout) :: ifile
      integer, intent(in) :: unit
      integer, intent(out), optional :: iostat
    end subroutine ifile_append_from_unit
    module subroutine ifile_write (ifile, unit, iostat)
      type(ifile_t), intent(in) :: ifile
      integer, intent(in), optional :: unit
      integer, intent(out), optional :: iostat
    end subroutine ifile_write
    module subroutine ifile_to_string_array (ifile, string)
      type(ifile_t), intent(in) :: ifile
      type(string_t), dimension(:), intent(inout), allocatable :: string
    end subroutine ifile_to_string_array
    module function ifile_get_length (ifile) result (length)
      integer :: length
      type(ifile_t), intent(in) :: ifile
    end function ifile_get_length
    module subroutine line_init (line, ifile, back)
      type(line_p), intent(inout) :: line
      type(ifile_t), intent(in) :: ifile
      logical, intent(in), optional :: back
    end subroutine line_init
    module subroutine line_final (line)
      type(line_p), intent(inout) :: line
    end subroutine line_final
    module subroutine line_advance (line)
      type(line_p), intent(inout) :: line
    end subroutine line_advance
    module subroutine line_backspace (line)
      type(line_p), intent(inout) :: line
    end subroutine line_backspace
    module function line_is_associated (line) result (ok)
      logical :: ok
      type(line_p), intent(in) :: line
    end function line_is_associated
    module function line_get_string (line) result (string)
      type(string_t) :: string
      type(line_p), intent(in) :: line
    end function line_get_string
    module function line_get_string_advance (line) result (string)
      type(string_t) :: string
      type(line_p), intent(inout) :: line
    end function line_get_string_advance
    module function line_get_index (line) result (index)
      integer :: index
      type(line_p), intent(in) :: line
    end function line_get_index
    module function line_get_length (line) result (length)
      integer :: length
      type(line_p), intent(in) :: line
    end function line_get_length
  end interface

end module ifiles
