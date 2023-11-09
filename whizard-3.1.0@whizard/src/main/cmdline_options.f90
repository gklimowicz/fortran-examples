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

module cmdline_options

  use iso_varying_string, string_t => varying_string
  use diagnostics

  implicit none
  private

  public :: init_options
  public :: no_option_value
  public :: get_option_value

  integer, parameter :: CMDLINE_ARG_LEN = 1000

  abstract interface
     subroutine msg
     end subroutine msg
  end interface

  procedure (msg), pointer :: print_usage => null ()

contains

  subroutine init_options (usage_msg)
    procedure (msg) :: usage_msg
    print_usage => usage_msg
  end subroutine init_options

  subroutine no_option_value (option, value)
    type(string_t), intent(in) :: option, value
    if (value /= "") then
       call msg_error (" Option '" // char (option) // "' should have no value")
    end if
  end subroutine no_option_value

  function get_option_value (i, option, value) result (string)
    type(string_t) :: string
    integer, intent(inout) :: i
    type(string_t), intent(in) :: option
    type(string_t), intent(in), optional :: value
    character(CMDLINE_ARG_LEN) :: arg_value
    integer :: arg_len, arg_status
    logical :: has_value
    if (present (value)) then
       has_value = value /= ""
    else
       has_value = .false.
    end if
    if (has_value) then
       call unquote_value (i, option, value, string)
    else
       i = i + 1
       call get_command_argument (i, arg_value, arg_len, arg_status)
       select case (arg_status)
       case (0)
       case (-1)
          call msg_error (" Option value truncated: '" // arg_value // "'")
       case default
          call print_usage ()
          call msg_fatal (" Option '" // char (option) // "' needs a value")
       end select
       select case (arg_value(1:1))
       case ("-")
          call print_usage ()
          call msg_fatal (" Option '" // char (option) // "' needs a value")
       end select
       call unquote_value (i, option, var_str (trim (arg_value)), string)
    end if
  end function get_option_value

  subroutine unquote_value (i, option, value, string)
    integer, intent(inout) :: i
    type(string_t), intent(in) :: option
    type(string_t), intent(in) :: value
    type(string_t), intent(out) :: string
    character(1) :: quote
    character(CMDLINE_ARG_LEN) :: arg_value
    integer :: arg_len, arg_status
    quote = extract (value, 1, 1)
    select case (quote)
    case ("'", '"')
       string = ""
       arg_value = extract (value, 2)
       arg_len = len_trim (value)
       APPEND_QUOTED: do
          if (extract (arg_value, arg_len, arg_len) == quote) then
             string = string // " " // extract (arg_value, 1, arg_len-1)
             exit APPEND_QUOTED
          else
             string = string // " " // trim (arg_value)
             i = i + 1
             call get_command_argument (i, arg_value, arg_len, arg_status)
             select case (arg_status)
             case (0)
             case (-1)
                call msg_error (" Quoted option value truncated: '" &
                     // char (string) // "'")
             case default
                call print_usage ()
                call msg_fatal (" Option '" // char (option) &
                     // "': unterminated quoted value")
             end select
          end if
       end do APPEND_QUOTED
    case default
       string = value
    end select
  end subroutine unquote_value

end module cmdline_options

