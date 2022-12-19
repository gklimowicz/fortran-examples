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

module user_files

  use iso_varying_string, string_t => varying_string
  use ifiles

  implicit none
  private

  public :: file_list_t
  public :: file_list_final
  public :: file_list_is_open
  public :: file_list_get_unit
  public :: file_list_open
  public :: file_list_close
  public :: file_list_write
  public :: file_list_write_analysis

  type :: file_t
     private
     type(string_t) :: name
     integer :: unit = -1
     logical :: reading = .false.
     logical :: writing = .false.
     type(file_t), pointer :: prev => null ()
     type(file_t), pointer :: next => null ()
  end type file_t

  type :: file_list_t
     type(file_t), pointer :: first => null ()
     type(file_t), pointer :: last => null ()
  end type file_list_t


  interface file_list_write
     module procedure file_list_write_string
     module procedure file_list_write_ifile
  end interface

  interface
    module subroutine file_list_final (file_list)
      type(file_list_t), intent(inout) :: file_list
    end subroutine file_list_final
    module function file_list_is_open (file_list, name, action) result (flag)
      logical :: flag
      type(file_list_t), intent(in) :: file_list
      type(string_t), intent(in) :: name
      character(len=*), intent(in) :: action
    end function file_list_is_open
    module function file_list_get_unit (file_list, name) result (unit)
      integer :: unit
      type(file_list_t), intent(in) :: file_list
      type(string_t), intent(in) :: name
    end function file_list_get_unit
    module subroutine file_list_open (file_list, name, action, status, position)
      type(file_list_t), intent(inout) :: file_list
      type(string_t), intent(in) :: name
      character(len=*), intent(in) :: action, status, position
    end subroutine file_list_open
    module subroutine file_list_close (file_list, name)
      type(file_list_t), intent(inout) :: file_list
      type(string_t), intent(in) :: name
    end subroutine file_list_close
    module subroutine file_list_write_ifile (file_list, name, ifile)
      type(file_list_t), intent(in) :: file_list
      type(string_t), intent(in) :: name
      type(ifile_t), intent(in) :: ifile
    end subroutine file_list_write_ifile
    module subroutine file_list_write_string &
         (file_list, name, string, advancing)
      type(file_list_t), intent(in) :: file_list
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: string
      logical, intent(in), optional :: advancing
    end subroutine file_list_write_string
    module subroutine file_list_write_analysis (file_list, name, tag)
      type(file_list_t), intent(in) :: file_list
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: tag
    end subroutine file_list_write_analysis
  end interface

end module user_files
