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

module prc_core_def

  use iso_varying_string, string_t => varying_string
  use process_constants
  use prclib_interfaces

  implicit none
  private

  public :: prc_core_def_t
  public :: prc_template_t
  public :: allocate_core_def
  public :: prc_core_driver_t
  public :: process_driver_internal_t

  type, abstract :: prc_core_def_t
     class(prc_writer_t), allocatable :: writer
   contains
     procedure (prc_core_def_get_string), nopass, deferred :: type_string
     procedure (prc_core_def_write), deferred :: write
     procedure (prc_core_def_read), deferred :: read
     procedure :: set_md5sum => prc_core_def_set_md5sum
     procedure(prc_core_def_allocate_driver), deferred :: allocate_driver
     procedure, nopass :: needs_code => prc_core_def_needs_code
     procedure(prc_core_def_get_features), nopass, deferred &
          :: get_features
     procedure(prc_core_def_connect), deferred :: connect
  end type prc_core_def_t

  type :: prc_template_t
     class(prc_core_def_t), allocatable :: core_def
  end type prc_template_t

  type, abstract :: prc_core_driver_t
   contains
     procedure(prc_core_driver_type_name), nopass, deferred :: type_name
  end type prc_core_driver_t

  type, extends (prc_core_driver_t), abstract :: process_driver_internal_t
   contains
     procedure(process_driver_fill_constants), deferred :: fill_constants
  end type process_driver_internal_t


  abstract interface
     function prc_core_def_get_string () result (string)
       import
       type(string_t) :: string
     end function prc_core_def_get_string
  end interface

  abstract interface
     subroutine prc_core_def_write (object, unit)
       import
       class(prc_core_def_t), intent(in) :: object
       integer, intent(in) :: unit
     end subroutine prc_core_def_write
  end interface

  abstract interface
     subroutine prc_core_def_read (object, unit)
       import
       class(prc_core_def_t), intent(out) :: object
       integer, intent(in) :: unit
     end subroutine prc_core_def_read
  end interface

  abstract interface
     subroutine prc_core_def_allocate_driver (object, driver, basename)
       import
       class(prc_core_def_t), intent(in) :: object
       class(prc_core_driver_t), intent(out), allocatable :: driver
       type(string_t), intent(in) :: basename
     end subroutine prc_core_def_allocate_driver
  end interface

  abstract interface
     subroutine prc_core_def_get_features (features)
       import
       type(string_t), dimension(:), allocatable, intent(out) :: features
     end subroutine prc_core_def_get_features
  end interface

  abstract interface
     subroutine prc_core_def_connect (def, lib_driver, i, proc_driver)
       import
       class(prc_core_def_t), intent(in) :: def
       class(prclib_driver_t), intent(in) :: lib_driver
       integer, intent(in) :: i
       class(prc_core_driver_t), intent(inout) :: proc_driver
     end subroutine prc_core_def_connect
  end interface

  abstract interface
     function prc_core_driver_type_name () result (type)
       import
       type(string_t) :: type
     end function prc_core_driver_type_name
  end interface

  abstract interface
     subroutine process_driver_fill_constants (driver, data)
       import
       class(process_driver_internal_t), intent(in) :: driver
       type(process_constants_t), intent(out) :: data
     end subroutine process_driver_fill_constants
  end interface


  interface
    module subroutine prc_core_def_set_md5sum (core_def, md5sum)
      class(prc_core_def_t), intent(inout) :: core_def
      character(32) :: md5sum
    end subroutine prc_core_def_set_md5sum
    module function prc_core_def_needs_code () result (flag)
      logical :: flag
    end function prc_core_def_needs_code
    module subroutine allocate_core_def (template, name, core_def)
      type(prc_template_t), dimension(:), intent(in) :: template
      type(string_t), intent(in) :: name
      class(prc_core_def_t), allocatable :: core_def
    end subroutine allocate_core_def
  end interface

end module prc_core_def
