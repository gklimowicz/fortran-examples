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

module formats

  use, intrinsic :: iso_c_binding

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: sprintf_arg_t
  public :: sprintf_arg_init
  public :: sprintf

  integer, parameter, public :: ARGTYPE_NONE = 0
  integer, parameter, public :: ARGTYPE_LOG = 1
  integer, parameter, public :: ARGTYPE_INT = 2
  integer, parameter, public :: ARGTYPE_REAL = 3
  integer, parameter, public :: ARGTYPE_STR = 4


  type :: sprintf_arg_t
    private
    integer :: type = ARGTYPE_NONE
    integer(c_int), dimension(:), allocatable :: ival
    real(c_double), dimension(:), allocatable :: rval
    character(c_char), dimension(:), allocatable :: sval
  end type sprintf_arg_t

  type :: sprintf_interface_t
    private
    character(c_char), dimension(:), allocatable :: input_fmt
    type(sprintf_arg_t) :: arg
    character(c_char), dimension(:), allocatable :: output_str
    integer :: output_str_len = 0
  end type sprintf_interface_t


  interface sprintf_arg_init
     module procedure sprintf_arg_init_log
     module procedure sprintf_arg_init_int
     module procedure sprintf_arg_init_real
     module procedure sprintf_arg_init_str
  end interface

  interface
    function sprintf_none (str, fmt) result (stat) bind(C)
      use iso_c_binding !NODEP!
      integer(c_int) :: stat
      character(c_char), dimension(*), intent(inout) :: str
      character(c_char), dimension(*), intent(in) :: fmt
    end function sprintf_none
  end interface

  interface
    function sprintf_int (str, fmt, val) result (stat) bind(C)
      use iso_c_binding !NODEP!
      integer(c_int) :: stat
      character(c_char), dimension(*), intent(inout) :: str
      character(c_char), dimension(*), intent(in) :: fmt
      integer(c_int), value :: val
    end function sprintf_int
  end interface

  interface
    function sprintf_double (str, fmt, val) result (stat) bind(C)
      use iso_c_binding !NODEP!
      integer(c_int) :: stat
      character(c_char), dimension(*), intent(inout) :: str
      character(c_char), dimension(*), intent(in) :: fmt
      real(c_double), value :: val
    end function sprintf_double
  end interface

  interface
    function sprintf_str(str, fmt, val) result (stat) bind(C)
      use iso_c_binding !NODEP!
      integer(c_int) :: stat
      character(c_char), dimension(*), intent(inout) :: str
      character(c_char), dimension(*), intent(in) :: fmt
      character(c_char), dimension(*), intent(in) :: val
    end function sprintf_str
  end interface



  interface
    module subroutine sprintf_arg_init_log (arg, lval)
      type(sprintf_arg_t), intent(out) :: arg
      logical, intent(in) :: lval
    end subroutine sprintf_arg_init_log
    module subroutine sprintf_arg_init_int (arg, ival)
      type(sprintf_arg_t), intent(out) :: arg
      integer, intent(in) :: ival
    end subroutine sprintf_arg_init_int
    module subroutine sprintf_arg_init_real (arg, rval)
      type(sprintf_arg_t), intent(out) :: arg
      real(default), intent(in) :: rval
    end subroutine sprintf_arg_init_real
    module subroutine sprintf_arg_init_str (arg, sval)
      type(sprintf_arg_t), intent(out) :: arg
      type(string_t), intent(in) :: sval
    end subroutine sprintf_arg_init_str
    module function sprintf (fmt, arg) result (string)
      type(string_t) :: string
      type(string_t), intent(in) :: fmt
      type(sprintf_arg_t), dimension(:), intent(in) :: arg
    end function sprintf
  end interface

end module formats
