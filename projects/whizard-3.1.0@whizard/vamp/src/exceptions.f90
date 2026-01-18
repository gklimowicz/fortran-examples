! exceptions.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module exceptions
  use kinds
  implicit none
  private
  public :: handle_exception
  public :: raise_exception, clear_exception, gather_exceptions
  
  integer, public, parameter :: &
       EXC_NONE = 0, &
       EXC_INFO = 1, &
       EXC_WARN = 2, &
       EXC_ERROR = 3, &
       EXC_FATAL = 4
  integer, private, parameter :: EXC_DEFAULT = EXC_ERROR
  integer, private, parameter :: NAME_LENGTH = 64
  type, public :: exception
     integer :: level = EXC_NONE
     character(len=NAME_LENGTH) :: message = ""
     character(len=NAME_LENGTH) :: origin = ""
  end type exception
contains
  subroutine handle_exception (exc)
    type(exception), intent(inout) :: exc
    character(len=10) :: name
    if (exc%level > 0) then
       select case (exc%level)
          case (EXC_NONE)
             name = "(none)"
          case (EXC_INFO)
             name = "info"
          case (EXC_WARN)
             name = "warning"
          case (EXC_ERROR)
             name = "error"
          case (EXC_FATAL)
             name = "fatal"
          case default
             name = "invalid"
       end select
       print *, trim (exc%origin), ": ", trim(name), ": ", trim (exc%message)
       if (exc%level >= EXC_FATAL) then
          print *, "terminated."
          stop
       end if
    end if
  end subroutine handle_exception
  elemental subroutine raise_exception (exc, level, origin, message)
    type(exception), intent(inout), optional :: exc
    integer, intent(in), optional :: level
    character(len=*), intent(in), optional :: origin, message
    integer :: local_level
    if (present (exc)) then
       if (present (level)) then
          local_level = level
       else
          local_level = EXC_DEFAULT
       end if
       if (exc%level < local_level) then
          exc%level = local_level
          if (present (origin)) then
             exc%origin = origin
          else
             exc%origin = "[vamp]"
          end if
          if (present (message)) then
             exc%message = message
          else
             exc%message = "[vamp]"
          end if
       end if
    end if
  end subroutine raise_exception
  elemental subroutine clear_exception (exc)
    type(exception), intent(inout) :: exc
    exc%level = 0
    exc%message = ""
    exc%origin = ""
  end subroutine clear_exception
  pure subroutine gather_exceptions (exc, excs)
    type(exception), intent(inout) :: exc
    type(exception), dimension(:), intent(in) :: excs
    integer :: i
    i = sum (maxloc (excs%level))
    if (exc%level < excs(i)%level) then
       call raise_exception (exc, excs(i)%level, excs(i)%origin, &
                             excs(i)%message)
    end if
  end subroutine gather_exceptions
end module exceptions
