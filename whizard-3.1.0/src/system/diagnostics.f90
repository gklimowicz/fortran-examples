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

module diagnostics

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  use system_defs, only: BUFFER_SIZE, MAX_ERRORS

  implicit none
  private

  public :: RESULT, DEBUG, DEBUG2
  public :: d_area
  public :: D_PARTICLES, D_EVENTS, D_SHOWER, D_MODEL_F, &
       D_MATCHING, D_TRANSFORMS, D_SUBTRACTION, D_VIRTUAL, D_THRESHOLD, &
       D_PHASESPACE, D_MISMATCH, D_ME_METHODS, D_PROCESS_INTEGRATION, &
       D_TAUOLA, D_CORE, D_VAMP2, D_MPI, D_QFT, D_BEAMS, D_REAL, D_FLAVOR
  public :: msg_level
  public :: set_debug_levels
  public :: set_debug2_levels
  public :: term_col
  public :: mask_fatal_errors
  public :: handle_fatal_errors
  public :: msg_count
  public :: msg_list_clear
  public :: msg_summary
  public :: msg_listing
  public :: msg_buffer
  public :: create_col_string
  public :: msg_terminate
  public :: msg_bug, msg_fatal, msg_error, msg_warning
  public :: msg_message, msg_result
  public :: msg_debug
  public :: msg_print_color
  public :: msg_debug2
  public :: debug_active
  public :: debug2_active
  public :: msg_show_progress
  public :: exit
  public :: msg_banner
  public :: logging
  public :: logfile_init
  public :: logfile_final
  public :: logfile_unit
  public :: expect_record
  public :: expect_clear
  public :: expect_summary
  public :: int2string
  public :: int2char
  public :: int2fixed
  public :: real2string
  public :: real2char
  public :: real2fixed
  public :: wo_sigint
  public :: wo_sigterm
  public :: wo_sigxcpu
  public :: wo_sigxfsz

  public :: mask_term_signals
  public :: release_term_signals
  public :: signal_is_pending
  public :: terminate_now_if_signal
  public :: single_event
  public :: terminate_now_if_single_event

  integer, parameter :: TERMINATE=-2, BUG=-1, FATAL=1, &
       ERROR=2, WARNING=3, MESSAGE=4, RESULT=5, &
       DEBUG=6, DEBUG2=7
  integer, parameter :: D_ALL=0, D_PARTICLES=1, D_EVENTS=2, &
       D_SHOWER=3, D_MODEL_F=4, &
       D_MATCHING=5, D_TRANSFORMS=6, &
       D_SUBTRACTION=7, D_VIRTUAL=8, D_THRESHOLD=9, D_PHASESPACE=10, &
       D_MISMATCH=11, D_ME_METHODS=12, D_PROCESS_INTEGRATION=13, &
       D_TAUOLA=14, D_CORE=15, D_VAMP2 = 16, D_MPI = 17, D_QFT = 18, &
       D_BEAMS=19, D_REAL=20, D_FLAVOR=21, D_LAST=21
  integer, parameter, public :: COL_UNDEFINED = -1
  integer, parameter, public :: COL_GREY = 90, COL_PEACH = 91, COL_LIGHT_GREEN = 92, &
     COL_LIGHT_YELLOW = 93, COL_LIGHT_BLUE = 94, COL_PINK = 95, &
     COL_LIGHT_AQUA = 96, COL_PEARL_WHITE = 97, COL_BLACK = 30, &
     COL_RED = 31, COL_GREEN = 32, COL_YELLOW = 33, COL_BLUE = 34, &
     COL_PURPLE = 35, COL_AQUA = 36

  integer, parameter, public :: TERM_STOP = 0, TERM_EXIT = 1, TERM_CRASH = 2

  type :: terminal_color_t
     integer :: color = COL_UNDEFINED
  contains
  
  end type terminal_color_t

  type :: string_list
     character(len=BUFFER_SIZE) :: string
     type(string_list), pointer :: next
  end type string_list
  type :: string_list_pointer
     type(string_list), pointer :: first, last
  end type string_list_pointer


  integer, save, dimension(D_ALL:D_LAST) :: msg_level = RESULT
  logical, save :: mask_fatal_errors = .false.
  integer, save :: handle_fatal_errors = TERM_EXIT
  integer, dimension(TERMINATE:WARNING), save :: msg_count = 0
  type(string_list_pointer), dimension(TERMINATE:WARNING), save :: &
       & msg_list = string_list_pointer (null(), null())
  character(len=BUFFER_SIZE), save :: msg_buffer = " "
  integer, save :: log_unit = -1
  logical, target, save :: logging = .false.
  integer, save :: expect_total = 0
  integer, save :: expect_failures = 0

  integer(c_int), bind(C), volatile :: wo_sigint = 0
  integer(c_int), bind(C), volatile :: wo_sigterm = 0
  integer(c_int), bind(C), volatile :: wo_sigxcpu = 0
  integer(c_int), bind(C), volatile :: wo_sigxfsz = 0

  logical :: single_event = .false.

  interface d_area
     module procedure d_area_of_string
     module procedure d_area_to_string
  end interface
  interface term_col
     module procedure term_col_int
     module procedure term_col_char
  end interface term_col

  interface msg_debug
     module procedure msg_debug_none
     module procedure msg_debug_logical
     module procedure msg_debug_integer
     module procedure msg_debug_real
     module procedure msg_debug_complex
     module procedure msg_debug_string
  end interface
  interface msg_print_color
     module procedure msg_print_color_none
     module procedure msg_print_color_logical
     module procedure msg_print_color_integer
     module procedure msg_print_color_real
  end interface
  interface msg_debug2
     module procedure msg_debug2_none
     module procedure msg_debug2_logical
     module procedure msg_debug2_integer
     module procedure msg_debug2_real
     module procedure msg_debug2_complex
     module procedure msg_debug2_string
  end interface
  interface
     subroutine exit (status) bind (C)
       use iso_c_binding !NODEP!
       integer(c_int), value :: status
     end subroutine exit
  end interface

  interface real2string
     module procedure real2string_list, real2string_fmt
  end interface
  interface real2char
     module procedure real2char_list, real2char_fmt
  end interface
  interface
     integer(c_int) function wo_mask_sigint () bind(C)
       import
     end function wo_mask_sigint
  end interface
  interface
     integer(c_int) function wo_mask_sigterm () bind(C)
       import
     end function wo_mask_sigterm
  end interface
  interface
     integer(c_int) function wo_mask_sigxcpu () bind(C)
       import
     end function wo_mask_sigxcpu
  end interface
  interface
     integer(c_int) function wo_mask_sigxfsz () bind(C)
       import
     end function wo_mask_sigxfsz
  end interface

  interface
     integer(c_int) function wo_release_sigint () bind(C)
       import
     end function wo_release_sigint
  end interface
  interface
     integer(c_int) function wo_release_sigterm () bind(C)
       import
     end function wo_release_sigterm
  end interface
  interface
     integer(c_int) function wo_release_sigxcpu () bind(C)
       import
     end function wo_release_sigxcpu
  end interface
  interface
     integer(c_int) function wo_release_sigxfsz () bind(C)
       import
     end function wo_release_sigxfsz
  end interface


  interface
    module function d_area_of_string (string) result (i)
      integer :: i
      type(string_t), intent(in) :: string
    end function d_area_of_string

    elemental module function d_area_to_string (i) result (string)
      type(string_t) :: string
      integer, intent(in) :: i
    end function d_area_to_string
    module subroutine set_debug_levels (area_str)
      type(string_t), intent(in) :: area_str
    end subroutine set_debug_levels
    module subroutine set_debug2_levels (area_str)
      type(string_t), intent(in) :: area_str
    end subroutine set_debug2_levels
    module function term_col_int (col_int) result (color)
      type(terminal_color_t) :: color
      integer, intent(in) :: col_int
    end function term_col_int

    module function term_col_char (col_char) result (color)
      type(terminal_color_t) :: color
      character(len=*), intent(in) :: col_char
    end function term_col_char
    module subroutine msg_list_clear
    end subroutine msg_list_clear
    module subroutine msg_summary (unit)
      integer, intent(in), optional :: unit
    end subroutine msg_summary
    module subroutine msg_listing (level, unit, prefix)
      integer, intent(in) :: level
      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: prefix
    end subroutine msg_listing
    module function create_col_string (color) result (col_string)
       type(string_t) :: col_string
       integer, intent(in) :: color
    end function create_col_string
  module subroutine msg_terminate (string, unit, quit_code)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    integer, intent(in), optional :: quit_code
  end subroutine msg_terminate

  module subroutine msg_bug (string, arr, unit)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
  end subroutine msg_bug

  recursive module subroutine msg_fatal (string, arr, unit)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
  end subroutine msg_fatal

  module subroutine msg_error (string, arr, unit)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
  end subroutine msg_error

  module subroutine msg_warning (string, arr, unit, color)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
    type(terminal_color_t), intent(in), optional :: color
  end subroutine msg_warning

  module subroutine msg_message (string, unit, arr, logfile, color)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
    logical, intent(in), optional :: logfile
    type(terminal_color_t), intent(in), optional :: color
  end subroutine msg_message

  module subroutine msg_result (string, arr, unit, logfile, color)
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: string
    type(string_t), dimension(:), intent(in), optional :: arr
    logical, intent(in), optional :: logfile
    type(terminal_color_t), intent(in), optional :: color
  end subroutine msg_result
    module subroutine msg_debug_none (area, string, color)
      integer, intent(in) :: area
      character(len=*), intent(in), optional :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug_none

    module subroutine msg_debug_logical (area, string, value, color)
      logical, intent(in) :: value
      integer, intent(in) :: area
      character(len=*), intent(in) :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug_logical

    module subroutine msg_debug_integer (area, string, value, color)
      integer, intent(in) :: value
      integer, intent(in) :: area
      character(len=*), intent(in) :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug_integer

    module subroutine msg_debug_real (area, string, value, color)
      real(default), intent(in) :: value
      integer, intent(in) :: area
      character(len=*), intent(in) :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug_real

    module subroutine msg_debug_complex (area, string, value, color)
      complex(default), intent(in) :: value
      integer, intent(in) :: area
      character(len=*), intent(in) :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug_complex

    module subroutine msg_debug_string (area, string, value, color)
     type(string_t), intent(in) :: value
      integer, intent(in) :: area
      character(len=*), intent(in) :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug_string
    module subroutine msg_print_color_none (string, color)
      character(len=*), intent(in) :: string
      !!!type(terminal_color_t), intent(in) :: color
      integer, intent(in) :: color
    end subroutine msg_print_color_none

    module subroutine msg_print_color_logical (string, value, color)
      character(len=*), intent(in) :: string
      logical, intent(in) :: value
      integer, intent(in) :: color
    end subroutine msg_print_color_logical

    module subroutine msg_print_color_integer (string, value, color)
      character(len=*), intent(in) :: string
      integer, intent(in) :: value
      integer, intent(in) :: color
    end subroutine msg_print_color_integer

    module subroutine msg_print_color_real (string, value, color)
      character(len=*), intent(in) :: string
      real(default), intent(in) :: value
      integer, intent(in) :: color
    end subroutine msg_print_color_real
    module subroutine msg_debug2_none (area, string, color)
      integer, intent(in) :: area
      character(len=*), intent(in), optional :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug2_none

    module subroutine msg_debug2_logical (area, string, value, color)
      logical, intent(in) :: value
      integer, intent(in) :: area
      character(len=*), intent(in) :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug2_logical

    module subroutine msg_debug2_integer (area, string, value, color)
      integer, intent(in) :: value
      integer, intent(in) :: area
      character(len=*), intent(in) :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug2_integer

    module subroutine msg_debug2_real (area, string, value, color)
      real(default), intent(in) :: value
      integer, intent(in) :: area
      character(len=*), intent(in) :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug2_real

    module subroutine msg_debug2_complex (area, string, value, color)
      complex(default), intent(in) :: value
      integer, intent(in) :: area
      character(len=*), intent(in) :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug2_complex

    module subroutine msg_debug2_string (area, string, value, color)
      type(string_t), intent(in) :: value
      integer, intent(in) :: area
      character(len=*), intent(in) :: string
      type(terminal_color_t), intent(in), optional :: color
    end subroutine msg_debug2_string
    elemental module function debug_active (area) result (active)
      logical :: active
      integer, intent(in) :: area
    end function debug_active
    elemental module function debug2_active (area) result (active)
      logical :: active
      integer, intent(in) :: area
    end function debug2_active
    module subroutine msg_show_progress (i_call, n_calls)
      integer, intent(in) :: i_call, n_calls
    end subroutine msg_show_progress
    module subroutine msg_banner (unit)
      integer, intent(in), optional :: unit
    end subroutine msg_banner
    module subroutine logfile_init (filename)
      type(string_t), intent(in) :: filename
    end subroutine logfile_init
    module subroutine logfile_final ()
    end subroutine logfile_final
    module function logfile_unit (unit, logfile)
      integer :: logfile_unit
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: logfile
    end function logfile_unit
    module subroutine expect_record (success)
      logical, intent(in) :: success
    end subroutine expect_record
    module subroutine expect_clear ()
    end subroutine expect_clear
    module subroutine expect_summary (unit, force)
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: force
    end subroutine expect_summary
    pure module function int2fixed (i) result (c)
      integer, intent(in) :: i
      character(200) :: c
    end function int2fixed

    pure module function int2string (i) result (s)
      integer, intent(in) :: i
      type (string_t) :: s
    end function int2string

    pure module function int2char (i) result (c)
      integer, intent(in) :: i
      character(len (trim (int2fixed (i)))) :: c
    end function int2char
    pure module function real2fixed (x, fmt) result (c)
      real(default), intent(in) :: x
      character(*), intent(in), optional :: fmt
      character(200) :: c
    end function real2fixed

    pure module function real2fixed_fmt (x, fmt) result (c)
      real(default), intent(in) :: x
      character(*), intent(in) :: fmt
      character(200) :: c
    end function real2fixed_fmt

    pure module function real2string_list (x) result (s)
      real(default), intent(in) :: x
      type(string_t) :: s
    end function real2string_list

    pure module function real2string_fmt (x, fmt) result (s)
      real(default), intent(in) :: x
      character(*), intent(in) :: fmt
      type(string_t) :: s
    end function real2string_fmt

    pure module function real2char_list (x) result (c)
      real(default), intent(in) :: x
      character(len_trim (real2fixed (x))) :: c
    end function real2char_list

    pure module function real2char_fmt (x, fmt) result (c)
      real(default), intent(in) :: x
      character(*), intent(in) :: fmt
      character(len_trim (real2fixed_fmt (x, fmt))) :: c
    end function real2char_fmt
    module subroutine mask_term_signals ()
    end subroutine mask_term_signals
    module subroutine release_term_signals ()
    end subroutine release_term_signals
    module function signal_is_pending () result (flag)
      logical :: flag
    end function signal_is_pending
    module subroutine terminate_now_if_signal ()
    end subroutine terminate_now_if_signal
    module subroutine terminate_now_if_single_event ()
    end subroutine terminate_now_if_single_event
  end interface

end module diagnostics

  subroutine fatal_force_crash ()
    use diagnostics, only: handle_fatal_errors, TERM_CRASH !NODEP!
    implicit none
    handle_fatal_errors = TERM_CRASH
  end subroutine fatal_force_crash

  subroutine fatal_force_exit ()
    use diagnostics, only: handle_fatal_errors, TERM_EXIT !NODEP!
    implicit none
    handle_fatal_errors = TERM_EXIT
  end subroutine fatal_force_exit

  subroutine fatal_force_stop ()
    use diagnostics, only: handle_fatal_errors, TERM_STOP !NODEP!
    implicit none
    handle_fatal_errors = TERM_STOP
  end subroutine fatal_force_stop

