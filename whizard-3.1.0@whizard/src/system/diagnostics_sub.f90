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

submodule (diagnostics) diagnostics_s

  use, intrinsic :: iso_fortran_env, only: output_unit !NODEP!

  use system_dependencies
  use debug_master, only: debug_on
  use string_utils, only: str
  use io_units

  implicit none

contains

  module function d_area_of_string (string) result (i)
    integer :: i
    type(string_t), intent(in) :: string
    select case (char (string))
    case ("particles")
       i = D_PARTICLES
    case ("events")
       i = D_EVENTS
    case ("shower")
       i = D_SHOWER
    case ("model_features")
       i = D_MODEL_F
    case ("matching")
       i = D_MATCHING
    case ("transforms")
       i = D_TRANSFORMS
    case ("subtraction")
       i = D_SUBTRACTION
    case ("virtual")
       i = D_VIRTUAL
    case ("threshold")
       i = D_THRESHOLD
    case ("phasespace")
       i = D_PHASESPACE
    case ("mismatch")
       i = D_MISMATCH
    case ("me_methods")
       i = D_ME_METHODS
    case ("process_integration")
       i = D_PROCESS_INTEGRATION
    case ("tauola")
       i = D_TAUOLA
    case ("core")
       i = D_CORE
    case ("vamp2")
       i = D_VAMP2
    case ("mpi")
       i = D_MPI
    case ("qft")
       i = D_QFT
    case ("beams")
       i = D_BEAMS
    case ("real")
       i = D_REAL
    case ("flavor")
       i = D_FLAVOR
    case ("all")
       i = D_ALL
    case default
       print "(A)", "Possible values for --debug are:"
       do i = 0, D_LAST
          print "(A)", char ('  ' // d_area_to_string(i))
       end do
       call msg_fatal ("Please use one of the listed areas")
    end select
  end function d_area_of_string

  elemental module function d_area_to_string (i) result (string)
    type(string_t) :: string
    integer, intent(in) :: i
    select case (i)
    case (D_PARTICLES)
       string = "particles"
    case (D_EVENTS)
       string = "events"
    case (D_SHOWER)
       string = "shower"
    case (D_MODEL_F)
       string = "model_features"
    case (D_MATCHING)
       string = "matching"
    case (D_TRANSFORMS)
       string = "transforms"
    case (D_SUBTRACTION)
       string = "subtraction"
    case (D_VIRTUAL)
       string = "virtual"
    case (D_THRESHOLD)
       string = "threshold"
    case (D_PHASESPACE)
       string = "phasespace"
    case (D_MISMATCH)
       string = "mismatch"
    case (D_ME_METHODS)
       string = "me_methods"
    case (D_PROCESS_INTEGRATION)
       string = "process_integration"
    case (D_TAUOLA)
       string = "tauola"
    case (D_CORE)
       string = "core"
    case (D_VAMP2)
       string = "vamp2"
    case (D_MPI)
       string = "mpi"
    case (D_QFT)
       string = "qft"
    case (D_BEAMS)
       string = "beams"
    case (D_REAL)
       string = "real"
    case (D_FLAVOR)
       string = "flavor"
    case (D_ALL)
       string = "all"
    case default
       string = "undefined"
    end select
  end function d_area_to_string

  module subroutine set_debug_levels (area_str)
    type(string_t), intent(in) :: area_str
    integer :: area
    if (.not. debug_on)  call msg_fatal ("Debugging options &
         &can be used only if configured with --enable-fc-debug")
    area = d_area (area_str)
    if (area == D_ALL) then
       msg_level = DEBUG
    else
       msg_level(area) = DEBUG
    end if
  end subroutine set_debug_levels

  module subroutine set_debug2_levels (area_str)
    type(string_t), intent(in) :: area_str
    integer :: area
    if (.not. debug_on)  call msg_fatal ("Debugging options &
         &can be used only if configured with --enable-fc-debug")
    area = d_area (area_str)
    if (area == D_ALL) then
       msg_level = DEBUG2
    else
       msg_level(area) = DEBUG2
    end if
  end subroutine set_debug2_levels

  module function term_col_int (col_int) result (color)
    type(terminal_color_t) :: color
    integer, intent(in) :: col_int
    color%color = col_int
  end function term_col_int

  module function term_col_char (col_char) result (color)
    type(terminal_color_t) :: color
    character(len=*), intent(in) :: col_char
    type(string_t) :: buf
    select case (col_char)
    case ('Grey')
       color%color = COL_GREY
    case ('Peach')
       color%color = COL_PEACH
    case ('Light Green')
       color%color = COL_LIGHT_GREEN
    case ('Light Yellow')
       color%color = COL_LIGHT_YELLOW
    case ('Light Blue')
       color%color = COL_LIGHT_BLUE
    case ('Pink')
       color%color = COL_PINK
    case ('Light Aqua')
       color%color = COL_LIGHT_AQUA
    case ('Pearl White')
       color%color = COL_PEARL_WHITE
    case ('Black')
       color%color = COL_BLACK
    case ('Red')
       color%color = COL_RED
    case ('Green')
       color%color = COL_GREEN
    case ('Yellow')
       color%color = COL_YELLOW
    case ('Blue')
       color%color = COL_BLUE
    case ('Purple')
       color%color = COL_PURPLE
    case ('Aqua')
       color%color = COL_AQUA
    case default
       buf = var_str ('Color ') // var_str (col_char) // var_str (' is not defined')
       call msg_warning (char (buf))
       color%color = COL_UNDEFINED
    end select
  end function term_col_char

  subroutine msg_add (level)
    integer, intent(in) :: level
    type(string_list), pointer :: message
    select case (level)
    case (TERMINATE:WARNING)
       allocate (message)
       message%string = msg_buffer
       nullify (message%next)
       if (.not.associated (msg_list(level)%first)) &
            & msg_list(level)%first => message
       if (associated (msg_list(level)%last)) &
            & msg_list(level)%last%next => message
       msg_list(level)%last => message
       msg_count(level) = msg_count(level) + 1
    end select
  end subroutine msg_add

  module subroutine msg_list_clear
    integer :: level
    type(string_list), pointer :: message
    do level = TERMINATE, WARNING
       do while (associated (msg_list(level)%first))
          message => msg_list(level)%first
          msg_list(level)%first => message%next
          deallocate (message)
       end do
       nullify (msg_list(level)%last)
    end do
    msg_count = 0
  end subroutine msg_list_clear

  module subroutine msg_summary (unit)
    integer, intent(in), optional :: unit
    call expect_summary (unit)
1   format (A,1x,I2,1x,A,I2,1x,A)
    if (msg_count(ERROR) > 0 .and. msg_count(WARNING) > 0) then
       write (msg_buffer, 1) "There were", &
            & msg_count(ERROR), "error(s) and  ", &
            & msg_count(WARNING), "warning(s)."
       call msg_message (unit=unit)
    else if (msg_count(ERROR) > 0) then
       write (msg_buffer, 1) "There were", &
            & msg_count(ERROR), "error(s) and no warnings."
       call msg_message (unit=unit)
    else if (msg_count(WARNING) > 0) then
       write (msg_buffer, 1) "There were no errors and  ", &
            & msg_count(WARNING), "warning(s)."
       call msg_message (unit=unit)
    end if
  end subroutine msg_summary

  module subroutine msg_listing (level, unit, prefix)
    integer, intent(in) :: level
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: prefix
    type(string_list), pointer :: message
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    if (present (unit))  u = unit
    message => msg_list(level)%first
    do while (associated (message))
       if (present (prefix)) then
          write (u, "(A)") prefix // trim (message%string)
       else
          write (u, "(A)") trim (message%string)
       end if
       message => message%next
    end do
    flush (u)
  end subroutine msg_listing

  subroutine buffer_clear
    msg_buffer = " "
  end subroutine buffer_clear

  module function create_col_string (color) result (col_string)
     type(string_t) :: col_string
     integer, intent(in) :: color
     character(2) :: buf
     write (buf, '(I2)') color
     col_string = var_str ("[") // var_str (buf) // var_str ("m")
  end function create_col_string

  subroutine message_print (level, string, str_arr, unit, logfile, area, color)
    integer, intent(in) :: level
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: str_arr
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: logfile
    integer, intent(in), optional :: area
    integer, intent(in), optional :: color
    type(string_t) :: col_string, prep_string, aux_string, head_footer, app_string
    integer :: lu, i, ar
    logical :: severe, is_error
    ar = D_ALL; if (present (area))  ar = area
    severe = .false.
    head_footer  = "******************************************************************************"
    aux_string = ""
    is_error = .false.
    app_string = ""
    select case (level)
    case (TERMINATE)
       prep_string = ""
    case (BUG)
       prep_string   = "*** WHIZARD BUG: "
       aux_string    = "***              "
       severe = .true.
       is_error = .true.
    case (FATAL)
       prep_string   = "*** FATAL ERROR: "
       aux_string    = "***              "
       severe = .true.
       is_error = .true.
    case (ERROR)
       prep_string = "*** ERROR: "
       aux_string  = "***        "
       is_error = .true.
    case (WARNING)
       prep_string = "Warning: "
    case (MESSAGE)
       prep_string = "| "
    case (DEBUG, DEBUG2)
       prep_string = "D: "
    case default
       prep_string = ""
    end select
    if (present (color)) then
       if (color > COL_UNDEFINED) then
          col_string = create_col_string (color)
          prep_string = achar(27) // col_string // prep_string
          app_string = app_string // achar(27) // "[0m"
       end if
    end if
    if (present(string))  msg_buffer = string
    lu = log_unit
    if (present(unit)) then
       if (unit /= output_unit) then
          if (severe) write (unit, "(A)") char(head_footer)
          if (is_error) write (unit, "(A)") char(head_footer)
          write (unit, "(A,A,A)") char(prep_string), trim(msg_buffer), &
               char(app_string)
          if (present (str_arr)) then
             do i = 1, size(str_arr)
                write (unit, "(A,A)") char(aux_string), char(trim(str_arr(i)))
             end do
          end if
          if (is_error) write (unit, "(A)") char(head_footer)
          if (severe) write (unit, "(A)") char(head_footer)
          flush (unit)
          lu = -1
       else if (level <= msg_level(ar)) then
          if (severe) print "(A)", char(head_footer)
          if (is_error) print "(A)", char(head_footer)
          print "(A,A,A)", char(prep_string), trim(msg_buffer), &
               char(app_string)
          if (present (str_arr)) then
             do i = 1, size(str_arr)
                print "(A,A)", char(aux_string), char(trim(str_arr(i)))
             end do
          end if
          if (is_error) print "(A)", char(head_footer)
          if (severe) print "(A)", char(head_footer)
          flush (output_unit)
          if (unit == log_unit)  lu = -1
       end if
    else if (level <= msg_level(ar)) then
       if (severe) print "(A)", char(head_footer)
       if (is_error) print "(A)", char(head_footer)
       print "(A,A,A)", char(prep_string), trim(msg_buffer), &
               char(app_string)
          if (present (str_arr)) then
             do i = 1, size(str_arr)
                print "(A,A)", char(aux_string), char(trim(str_arr(i)))
             end do
          end if
       if (is_error) print "(A)", char(head_footer)
       if (severe) print "(A)", char(head_footer)
       flush (output_unit)
    end if
    if (present (logfile)) then
       if (.not. logfile)  lu = -1
    end if
    if (logging .and. lu >= 0) then
       if (severe) write (lu, "(A)") char(head_footer)
       if (is_error) write (lu, "(A)") char(head_footer)
       write (lu, "(A,A,A)")  char(prep_string), trim(msg_buffer), &
               char(app_string)
       if (present (str_arr)) then
          do i = 1, size(str_arr)
             write (lu, "(A,A)") char(aux_string), char(trim(str_arr(i)))
          end do
       end if
       if (is_error) write (lu, "(A)") char(head_footer)
       if (severe) write (lu, "(A)") char(head_footer)
       flush (lu)
    end if
    call msg_add (level)
    call buffer_clear
  end subroutine message_print

  module subroutine msg_terminate (string, unit, quit_code)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    integer, intent(in), optional :: quit_code
    integer(c_int) :: return_code
    call release_term_signals ()
    if (present (quit_code)) then
       return_code = quit_code
    else
       return_code = 0
    end if
    if (present (string)) &
         call message_print (MESSAGE, string, unit=unit)
    call msg_summary (unit)
    if (return_code == 0 .and. expect_failures /= 0) then
       return_code = 5
       call message_print (MESSAGE, &
            "WHIZARD run finished with 'expect' failure(s).", unit=unit)
    else if (return_code == 7) then
       call message_print (MESSAGE, &
            "WHIZARD run finished with failed self-test.", unit=unit)
    else
       call message_print (MESSAGE, "WHIZARD run finished.", unit=unit)
    end if
    call message_print (0, &
         "|=============================================================================|", unit=unit)
    call logfile_final ()
    call msg_list_clear ()
    if (return_code /= 0) then
       call exit (return_code)
    else
       !!! Should implement WHIZARD exit code (currently only via C)
       call exit (0)
    end if
  end subroutine msg_terminate

  module subroutine msg_bug (string, arr, unit)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
    logical, pointer :: crash_ptr
    call message_print (BUG, string, arr, unit)
    call msg_summary (unit)
    select case (handle_fatal_errors)
    case (TERM_EXIT)
       call message_print (TERMINATE, "WHIZARD run aborted.", unit=unit)
       call exit (-1_c_int)
    case (TERM_CRASH)
       print *, "*** Intentional crash ***"
       crash_ptr => null ()
       print *, crash_ptr
    end select
    stop "WHIZARD run aborted."
  end subroutine msg_bug

  recursive module subroutine msg_fatal (string, arr, unit)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
    logical, pointer :: crash_ptr
    if (mask_fatal_errors) then
       call msg_error (string, arr, unit)
    else
       call message_print (FATAL, string, arr, unit)
       call msg_summary (unit)
       select case (handle_fatal_errors)
       case (TERM_EXIT)
          call message_print (TERMINATE, "WHIZARD run aborted.", unit=unit)
          call exit (1_c_int)
       case (TERM_CRASH)
          print *, "*** Intentional crash ***"
          crash_ptr => null ()
          print *, crash_ptr
       end select
       stop "WHIZARD run aborted."
    end if
  end subroutine msg_fatal

  module subroutine msg_error (string, arr, unit)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
    call message_print (ERROR, string, arr, unit)
    if (msg_count(ERROR) >= MAX_ERRORS) then
       mask_fatal_errors = .false.
       call msg_fatal (" Too many errors encountered.")
    else if (.not.present(unit) .and. .not.mask_fatal_errors)  then
       call message_print (MESSAGE, "            (WHIZARD run continues)")
    end if
  end subroutine msg_error

  module subroutine msg_warning (string, arr, unit, color)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
    type(terminal_color_t), intent(in), optional :: color
    integer :: cl
    cl = COL_UNDEFINED; if (present (color)) cl = color%color
    call message_print (level = WARNING, string = string, &
       str_arr = arr, unit = unit, color = cl)
  end subroutine msg_warning

  module subroutine msg_message (string, unit, arr, logfile, color)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
    logical, intent(in), optional :: logfile
    type(terminal_color_t), intent(in), optional :: color
    integer :: cl
    cl = COL_UNDEFINED; if (present (color)) cl = color%color
    call message_print (level = MESSAGE, &
       string = string, str_arr = arr, unit = unit, &
       logfile = logfile, color = cl)
  end subroutine msg_message

  module subroutine msg_result (string, arr, unit, logfile, color)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
    logical, intent(in), optional :: logfile
    type(terminal_color_t), intent(in), optional :: color
    integer :: cl
    cl = COL_UNDEFINED; if (present (color)) cl = color%color
    call message_print (level = RESULT, string = string, &
       str_arr = arr, unit = unit, logfile = logfile, color = cl)
  end subroutine msg_result

  module subroutine msg_debug_none (area, string, color)
    integer, intent(in) :: area
    character(len=*), intent(in), optional :: string
    type(terminal_color_t), intent(in), optional :: color
    integer :: cl
    if (debug_active (area)) then
       cl = COL_BLUE; if (present (color)) cl = color%color
       call message_print (DEBUG, string, unit = output_unit, &
            area = area, logfile = .false., color = cl)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug called with debug_on=.false.")
    end if
  end subroutine msg_debug_none

  module subroutine msg_debug_logical (area, string, value, color)
    logical, intent(in) :: value
    integer, intent(in) :: area
    character(len=*), intent(in) :: string
    type(terminal_color_t), intent(in), optional :: color
    character(len=64) :: buffer
    if (debug_active (area)) then
       write (buffer, *)  value
       call msg_debug_none (area, string // " = " // trim (buffer), &
            color = color)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug called with debug_on=.false.")
    end if
  end subroutine msg_debug_logical

  module subroutine msg_debug_integer (area, string, value, color)
    integer, intent(in) :: value
    integer, intent(in) :: area
    character(len=*), intent(in) :: string
    type(terminal_color_t), intent(in), optional :: color
    character(len=64) :: buffer
    if (debug_active (area)) then
       write (buffer, *)  value
       call msg_debug_none (area, string // " = " // trim (buffer), &
            color = color)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug called with debug_on=.false.")
    end if
  end subroutine msg_debug_integer

  module subroutine msg_debug_real (area, string, value, color)
    real(default), intent(in) :: value
    integer, intent(in) :: area
    character(len=*), intent(in) :: string
    type(terminal_color_t), intent(in), optional :: color
    character(len=64) :: buffer
    if (debug_active (area)) then
       write (buffer, *)  value
       call msg_debug_none (area, string // " = " // trim (buffer), &
            color = color)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug called with debug_on=.false.")
    end if
  end subroutine msg_debug_real

  module subroutine msg_debug_complex (area, string, value, color)
    complex(default), intent(in) :: value
    integer, intent(in) :: area
    character(len=*), intent(in) :: string
    type(terminal_color_t), intent(in), optional :: color
    character(len=64) :: buffer
    if (debug_active (area)) then
       write (buffer, *)  value
       call msg_debug_none (area, string // " = " // trim (buffer), &
            color = color)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug called with debug_on=.false.")
    end if
  end subroutine msg_debug_complex

  module subroutine msg_debug_string (area, string, value, color)
    type(string_t), intent(in) :: value
    integer, intent(in) :: area
    character(len=*), intent(in) :: string
    type(terminal_color_t), intent(in), optional :: color
    if (debug_active (area)) then
       call msg_debug_none (area, string // " = " // char (value), &
            color = color)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug called with debug_on=.false.")
    end if
  end subroutine msg_debug_string

  module subroutine msg_print_color_none (string, color)
    character(len=*), intent(in) :: string
    !!!type(terminal_color_t), intent(in) :: color
    integer, intent(in) :: color
    call message_print (0, string, color = color)
  end subroutine msg_print_color_none

  module subroutine msg_print_color_logical (string, value, color)
    character(len=*), intent(in) :: string
    logical, intent(in) :: value
    integer, intent(in) :: color
    call msg_print_color_none (char (string // " = " // str (value)), &
       color = color)
  end subroutine msg_print_color_logical

  module subroutine msg_print_color_integer (string, value, color)
    character(len=*), intent(in) :: string
    integer, intent(in) :: value
    integer, intent(in) :: color
    call msg_print_color_none (char (string // " = " // str (value)), &
       color = color)
  end subroutine msg_print_color_integer

  module subroutine msg_print_color_real (string, value, color)
    character(len=*), intent(in) :: string
    real(default), intent(in) :: value
    integer, intent(in) :: color
    call msg_print_color_none (char (string // " = " // str (value)), &
       color = color)
  end subroutine msg_print_color_real

  module subroutine msg_debug2_none (area, string, color)
    integer, intent(in) :: area
    character(len=*), intent(in), optional :: string
    type(terminal_color_t), intent(in), optional :: color
    integer :: cl
    if (debug2_active (area)) then
       cl = COL_BLUE; if (present (color)) cl = color%color
       call message_print (DEBUG2, string, unit = output_unit, &
            area = area, logfile = .false., color = cl)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug2 called with debug_on=.false.")
    end if
  end subroutine msg_debug2_none

  module subroutine msg_debug2_logical (area, string, value, color)
    logical, intent(in) :: value
    integer, intent(in) :: area
    character(len=*), intent(in) :: string
    type(terminal_color_t), intent(in), optional :: color
    character(len=64) :: buffer
    if (debug2_active (area)) then
       write (buffer, *)  value
       call msg_debug2_none (area, string // " = " // trim (buffer), &
            color = color)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug2 called with debug_on=.false.")
    end if
  end subroutine msg_debug2_logical

  module subroutine msg_debug2_integer (area, string, value, color)
    integer, intent(in) :: value
    integer, intent(in) :: area
    character(len=*), intent(in) :: string
    type(terminal_color_t), intent(in), optional :: color
    character(len=64) :: buffer
    if (debug2_active (area)) then
       write (buffer, *)  value
       call msg_debug2_none (area, string // " = " // trim (buffer), &
            color = color)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug2 called with debug_on=.false.")
    end if
  end subroutine msg_debug2_integer

  module subroutine msg_debug2_real (area, string, value, color)
    real(default), intent(in) :: value
    integer, intent(in) :: area
    character(len=*), intent(in) :: string
    type(terminal_color_t), intent(in), optional :: color
    character(len=64) :: buffer
    if (debug2_active (area)) then
       write (buffer, *)  value
       call msg_debug2_none (area, string // " = " // trim (buffer), &
            color = color)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug2 called with debug_on=.false.")
    end if
  end subroutine msg_debug2_real

  module subroutine msg_debug2_complex (area, string, value, color)
    complex(default), intent(in) :: value
    integer, intent(in) :: area
    character(len=*), intent(in) :: string
    type(terminal_color_t), intent(in), optional :: color
    character(len=64) :: buffer
    if (debug2_active (area)) then
       write (buffer, *)  value
       call msg_debug2_none (area, string // " = " // trim (buffer), &
            color = color)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug2 called with debug_on=.false.")
    end if
  end subroutine msg_debug2_complex

  module subroutine msg_debug2_string (area, string, value, color)
    type(string_t), intent(in) :: value
    integer, intent(in) :: area
    character(len=*), intent(in) :: string
    type(terminal_color_t), intent(in), optional :: color
    if (debug2_active (area)) then
       call msg_debug2_none (area, string // " = " // char (value), &
            color = color)
    else
       if (.not. debug_on)  call msg_bug ("msg_debug2 called with debug_on=.false.")
    end if
  end subroutine msg_debug2_string

  elemental module function debug_active (area) result (active)
    logical :: active
    integer, intent(in) :: area
    active = debug_on .and. msg_level(area) >= DEBUG
  end function debug_active

  elemental module function debug2_active (area) result (active)
    logical :: active
    integer, intent(in) :: area
    active = debug_on .and. msg_level(area) >= DEBUG2
  end function debug2_active

  module subroutine msg_show_progress (i_call, n_calls)
    integer, intent(in) :: i_call, n_calls
    real(default) :: progress
    integer, save :: next_check
    if (i_call == 1) next_check = 10
    progress = (i_call * 100._default) / n_calls
    if (progress >= next_check) then
       write (msg_buffer, "(F5.1,A)") progress, "%"
       call msg_message ()
       next_check = next_check + 10
    end if
  end subroutine msg_show_progress

  module subroutine msg_banner (unit)
    integer, intent(in), optional :: unit
    call message_print (0, "|=============================================================================|", unit=unit)
    call message_print (0, "|                                                                             |", unit=unit)
    call message_print (0, "|    WW             WW  WW   WW  WW  WWWWWW      WW      WWWWW    WWWW        |", unit=unit)
    call message_print (0, "|     WW    WW     WW   WW   WW  WW     WW      WWWW     WW  WW   WW  WW      |", unit=unit)
    call message_print (0, "|      WW  WW WW  WW    WWWWWWW  WW    WW      WW  WW    WWWWW    WW   WW     |", unit=unit)
    call message_print (0, "|       WWWW   WWWW     WW   WW  WW   WW      WWWWWWWW   WW  WW   WW  WW      |", unit=unit)
    call message_print (0, "|        WW     WW      WW   WW  WW  WWWWWW  WW      WW  WW   WW  WWWW        |", unit=unit)
    call message_print (0, "|                                                                             |", unit=unit)
    call message_print (0, "|                                                                             |", unit=unit)
    call message_print (0, "|                                        W                                    |", unit=unit)
    call message_print (0, "|                                       sW                                    |", unit=unit)
    call message_print (0, "|                                       WW                                    |", unit=unit)
    call message_print (0, "|                                      sWW                                    |", unit=unit)
    call message_print (0, "|                                      WWW                                    |", unit=unit)
    call message_print (0, "|                                     wWWW                                    |", unit=unit)
    call message_print (0, "|                                    wWWWW                                    |", unit=unit)
    call message_print (0, "|                                    WW WW                                    |", unit=unit)
    call message_print (0, "|                                    WW WW                                    |", unit=unit)
    call message_print (0, "|                                   wWW WW                                    |", unit=unit)
    call message_print (0, "|                                  wWW  WW                                    |", unit=unit)
    call message_print (0, "|                                  WW   WW                                    |", unit=unit)
    call message_print (0, "|                                  WW   WW                                    |", unit=unit)
    call message_print (0, "|                                 WW    WW                                    |", unit=unit)
    call message_print (0, "|                                 WW    WW                                    |", unit=unit)
    call message_print (0, "|                                WW     WW                                    |", unit=unit)
    call message_print (0, "|                                WW     WW                                    |", unit=unit)
    call message_print (0, "|           wwwwww              WW      WW                                    |", unit=unit)
    call message_print (0, "|              WWWWWww          WW      WW                                    |", unit=unit)
    call message_print (0, "|                 WWWWWwwwww   WW       WW                                    |", unit=unit)
    call message_print (0, "|                     wWWWwwwwwWW       WW                                    |", unit=unit)
    call message_print (0, "|                 wWWWWWWWWWWwWWW       WW                                    |", unit=unit)
    call message_print (0, "|                wWWWWW       wW        WWWWWWW                               |", unit=unit)
    call message_print (0, "|                  WWWW       wW        WW  wWWWWWWWwww                       |", unit=unit)
    call message_print (0, "|                   WWWW                      wWWWWWWWwwww                    |", unit=unit)
    call message_print (0, "|                     WWWW                      WWWW     WWw                  |", unit=unit)
    call message_print (0, "|                       WWWWww                   WWWW                         |", unit=unit)
    call message_print (0, "|                           WWWwwww              WWWW                         |", unit=unit)
    call message_print (0, "|                               wWWWWwww       wWWWWW                         |", unit=unit)
    call message_print (0, "|                                     WwwwwwwwwWWW                            |", unit=unit)
    call message_print (0, "|                                                                             |", unit=unit)
    call message_print (0, "|                                                                             |", unit=unit)
    call message_print (0, "|                                                                             |", unit=unit)
    call message_print (0, "|  by:   Wolfgang Kilian, Thorsten Ohl, Juergen Reuter                        |", unit=unit)
    call message_print (0, "|        with contributions from Christian Speckner                           |", unit=unit)
    call message_print (0, "|        Contact: <whizard@desy.de>                                           |", unit=unit)
    call message_print (0, "|                                                                             |", unit=unit)
    call message_print (0, "|  if you use WHIZARD please cite:                                            |", unit=unit)
    call message_print (0, "|        W. Kilian, T. Ohl, J. Reuter,  Eur.Phys.J.C71 (2011) 1742            |", unit=unit)
    call message_print (0, "|                                          [arXiv: 0708.4233 [hep-ph]]        |", unit=unit)
    call message_print (0, "|        M. Moretti, T. Ohl, J. Reuter, arXiv: hep-ph/0102195                 |", unit=unit)
    call message_print (0, "|                                                                             |", unit=unit)
    call message_print (0, "|=============================================================================|", unit=unit)
    call message_print (0, "|                               WHIZARD " // WHIZARD_VERSION, unit=unit)
    call message_print (0, "|=============================================================================|", unit=unit)
  end subroutine msg_banner

  module subroutine logfile_init (filename)
    type(string_t), intent(in) :: filename
    call msg_message ("Writing log to '" // char (filename) // "'")
    if (.not. logging)  call msg_message ("(Logging turned off.)")
    log_unit = free_unit ()
    open (file = char (filename), unit = log_unit, &
          action = "write", status = "replace")
  end subroutine logfile_init

  module subroutine logfile_final ()
    if (log_unit >= 0) then
       close (log_unit)
       log_unit = -1
    end if
  end subroutine logfile_final

  module function logfile_unit (unit, logfile)
    integer :: logfile_unit
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: logfile
    if (logging) then
       if (present (unit)) then
          if (unit == output_unit) then
             logfile_unit = log_unit
          else
             logfile_unit = -1
          end if
       else if (present (logfile)) then
          if (logfile) then
             logfile_unit = log_unit
          else
             logfile_unit = -1
          end if
       else
          logfile_unit = log_unit
       end if
    else
       logfile_unit = -1
    end if
  end function logfile_unit

  module subroutine expect_record (success)
    logical, intent(in) :: success
    expect_total = expect_total + 1
    if (.not. success)  expect_failures = expect_failures + 1
  end subroutine expect_record

  module subroutine expect_clear ()
    expect_total = 0
    expect_failures = 0
  end subroutine expect_clear

  module subroutine expect_summary (unit, force)
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: force
    logical :: force_output
    force_output = .false.;  if (present (force))  force_output = force
    if (expect_total /= 0 .or. force_output) then
       call msg_message ("Summary of value checks:", unit)
       write (msg_buffer, "(2x,A,1x,I0,1x,A,1x,A,1x,I0)") &
            "Failures:", expect_failures, "/", "Total:", expect_total
       call msg_message (unit=unit)
    end if
  end subroutine expect_summary

  pure module function int2fixed (i) result (c)
    integer, intent(in) :: i
    character(200) :: c
    c = ""
    write (c, *) i
    c = adjustl (c)
  end function int2fixed

  pure module function int2string (i) result (s)
    integer, intent(in) :: i
    type (string_t) :: s
    s = trim (int2fixed (i))
  end function int2string

  pure module function int2char (i) result (c)
    integer, intent(in) :: i
    character(len (trim (int2fixed (i)))) :: c
    c = int2fixed (i)
  end function int2char

  pure module function real2fixed (x, fmt) result (c)
    real(default), intent(in) :: x
    character(*), intent(in), optional :: fmt
    character(200) :: c
    c = ""
    write (c, *) x
    c = adjustl (c)
  end function real2fixed

  pure module function real2fixed_fmt (x, fmt) result (c)
    real(default), intent(in) :: x
    character(*), intent(in) :: fmt
    character(200) :: c
    c = ""
    write (c, fmt)  x
    c = adjustl (c)
  end function real2fixed_fmt

  pure module function real2string_list (x) result (s)
    real(default), intent(in) :: x
    type(string_t) :: s
    s = trim (real2fixed (x))
  end function real2string_list

  pure module function real2string_fmt (x, fmt) result (s)
    real(default), intent(in) :: x
    character(*), intent(in) :: fmt
    type(string_t) :: s
    s = trim (real2fixed_fmt (x, fmt))
  end function real2string_fmt

  pure module function real2char_list (x) result (c)
    real(default), intent(in) :: x
    character(len_trim (real2fixed (x))) :: c
    c = real2fixed (x)
  end function real2char_list

  pure module function real2char_fmt (x, fmt) result (c)
    real(default), intent(in) :: x
    character(*), intent(in) :: fmt
    character(len_trim (real2fixed_fmt (x, fmt))) :: c
    c = real2fixed_fmt (x, fmt)
  end function real2char_fmt

  module subroutine mask_term_signals ()
    logical :: ok
    wo_sigint = 0
    ok = wo_mask_sigint () == 0
    if (.not. ok)  call msg_error ("Masking SIGINT failed")
    wo_sigterm = 0
    ok = wo_mask_sigterm () == 0
    if (.not. ok)  call msg_error ("Masking SIGTERM failed")
    wo_sigxcpu = 0
    ok = wo_mask_sigxcpu () == 0
    if (.not. ok)  call msg_error ("Masking SIGXCPU failed")
    wo_sigxfsz = 0
    ok = wo_mask_sigxfsz () == 0
    if (.not. ok)  call msg_error ("Masking SIGXFSZ failed")
  end subroutine mask_term_signals

  module subroutine release_term_signals ()
    logical :: ok
    ok = wo_release_sigint () == 0
    if (.not. ok)  call msg_error ("Releasing SIGINT failed")
    ok = wo_release_sigterm () == 0
    if (.not. ok)  call msg_error ("Releasing SIGTERM failed")
    ok = wo_release_sigxcpu () == 0
    if (.not. ok)  call msg_error ("Releasing SIGXCPU failed")
    ok = wo_release_sigxfsz () == 0
    if (.not. ok)  call msg_error ("Releasing SIGXFSZ failed")
  end subroutine release_term_signals

  module function signal_is_pending () result (flag)
    logical :: flag
    flag = &
         wo_sigint /= 0 .or. &
         wo_sigterm /= 0 .or. &
         wo_sigxcpu /= 0 .or. &
         wo_sigxfsz /= 0
  end function signal_is_pending

  module subroutine terminate_now_if_signal ()
    if (wo_sigint /= 0) then
       call msg_terminate ("Signal SIGINT (keyboard interrupt) received.", &
          quit_code=int (wo_sigint))
    else if (wo_sigterm /= 0) then
       call msg_terminate ("Signal SIGTERM (termination signal) received.", &
          quit_code=int (wo_sigterm))
    else if (wo_sigxcpu /= 0) then
       call msg_terminate ("Signal SIGXCPU (CPU time limit exceeded) received.", &
          quit_code=int (wo_sigxcpu))
    else if (wo_sigxfsz /= 0) then
       call msg_terminate ("Signal SIGXFSZ (file size limit exceeded) received.", &
          quit_code=int (wo_sigxfsz))
    end if
  end subroutine terminate_now_if_signal

  module subroutine terminate_now_if_single_event ()
    integer, save :: n_calls = 0
    n_calls = n_calls + 1
    if (single_event .and. n_calls > 1) then
       call msg_terminate ("Stopping after one event", quit_code=0)
    end if
  end subroutine terminate_now_if_single_event


end submodule diagnostics_s

