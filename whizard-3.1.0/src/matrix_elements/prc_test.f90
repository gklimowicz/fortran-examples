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

module prc_test

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use particle_specifiers, only: new_prt_spec
  use process_constants
  use prclib_interfaces
  use prc_core_def
  use process_libraries

  implicit none
  private

  public :: prc_test_def_t
  public :: prc_test_t
  public :: prc_test_create_library

  type, extends (prc_core_def_t) :: prc_test_def_t
     type(string_t) :: model_name
     type(string_t), dimension(:), allocatable :: prt_in
     type(string_t), dimension(:), allocatable :: prt_out
   contains
     procedure, nopass :: type_string => prc_test_def_type_string
     procedure, nopass :: get_features => prc_test_def_get_features
     procedure :: init => prc_test_def_init
     procedure :: write => prc_test_def_write
     procedure :: read => prc_test_def_read
     procedure :: allocate_driver => prc_test_def_allocate_driver
     procedure :: connect => prc_test_def_connect
  end type prc_test_def_t

  type, extends (process_driver_internal_t) :: prc_test_t
     type(string_t) :: id
     type(string_t) :: model_name
     logical :: scattering = .true.
   contains
     procedure, nopass :: get_amplitude => prc_test_get_amplitude
     procedure, nopass :: type_name => prc_test_type_name
     procedure :: fill_constants => prc_test_fill_constants
  end type prc_test_t


  interface
    module function prc_test_def_type_string () result (string)
      type(string_t) :: string
    end function prc_test_def_type_string
    module subroutine prc_test_def_get_features (features)
      type(string_t), dimension(:), allocatable, intent(out) :: features
    end subroutine prc_test_def_get_features
    module subroutine prc_test_def_init (object, model_name, prt_in, prt_out)
      class(prc_test_def_t), intent(out) :: object
      type(string_t), intent(in) :: model_name
      type(string_t), dimension(:), intent(in) :: prt_in
      type(string_t), dimension(:), intent(in) :: prt_out
    end subroutine prc_test_def_init
    module subroutine prc_test_def_write (object, unit)
      class(prc_test_def_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine prc_test_def_write
    module subroutine prc_test_def_read (object, unit)
      class(prc_test_def_t), intent(out) :: object
      integer, intent(in) :: unit
    end subroutine prc_test_def_read
    module subroutine prc_test_def_connect (def, lib_driver, i, proc_driver)
      class(prc_test_def_t), intent(in) :: def
      class(prclib_driver_t), intent(in) :: lib_driver
      integer, intent(in) :: i
      class(prc_core_driver_t), intent(inout) :: proc_driver
    end subroutine prc_test_def_connect
    module function prc_test_get_amplitude (p) result (amp)
      complex(default) :: amp
      real(default), dimension(:,:), intent(in) :: p
    end function prc_test_get_amplitude
    module function prc_test_type_name () result (string)
      type(string_t) :: string
    end function prc_test_type_name
    module subroutine prc_test_fill_constants (driver, data)
      class(prc_test_t), intent(in) :: driver
      type(process_constants_t), intent(out) :: data
    end subroutine prc_test_fill_constants
  end interface

contains

  subroutine prc_test_def_allocate_driver (object, driver, basename)
    class(prc_test_def_t), intent(in) :: object
    class(prc_core_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    allocate (prc_test_t :: driver)
    select type (driver)
    type is (prc_test_t)
       driver%id = basename
       driver%model_name = object%model_name
       select case (size (object%prt_in))
       case (1);  driver%scattering = .false.
       case (2);  driver%scattering = .true.
       end select
    end select
  end subroutine prc_test_def_allocate_driver

  subroutine prc_test_create_library &
       (libname, lib, scattering, decay, procname1, procname2)
    type(string_t), intent(in) :: libname
    type(process_library_t), intent(out) :: lib
    logical, intent(in), optional :: scattering, decay
    type(string_t), intent(in), optional :: procname1, procname2
    type(string_t) :: model_name, procname
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry
    type(os_data_t) :: os_data
    logical :: sca, dec
    sca = .true.;   if (present (scattering))  sca = scattering
    dec = .false.;  if (present (decay))       dec = decay

    call os_data%init ()
    call lib%init (libname)
    model_name = "Test"

    if (sca) then
       if (present (procname1)) then
          procname = procname1
       else
          procname = libname
       end if
       allocate (prt_in (2), prt_out (2))
       prt_in  = [var_str ("s"), var_str ("s")]
       prt_out = [var_str ("s"), var_str ("s")]
       allocate (prc_test_def_t :: def)
       select type (def)
       type is (prc_test_def_t)
          call def%init (model_name, prt_in, prt_out)
       end select
       allocate (entry)
       call entry%init (procname, model_name = model_name, &
            n_in = 2, n_components = 1)
       call entry%import_component (1, n_out = size (prt_out), &
            prt_in  = new_prt_spec (prt_in), &
            prt_out = new_prt_spec (prt_out), &
            method  = var_str ("test_me"), &
            variant = def)
       call lib%append (entry)
    end if

    if (dec) then
       if (present (procname2)) then
          procname = procname2
       else
          procname = libname
       end if
       if (allocated (prt_in))  deallocate (prt_in, prt_out)
       allocate (prt_in (1), prt_out (2))
       prt_in  = [var_str ("s")]
       prt_out = [var_str ("f"), var_str ("fbar")]
       allocate (prc_test_def_t :: def)
       select type (def)
       type is (prc_test_def_t)
          call def%init (model_name, prt_in, prt_out)
       end select
       allocate (entry)
       call entry%init (procname, model_name = model_name, &
            n_in = 1, n_components = 1)
       call entry%import_component (1, n_out = size (prt_out), &
            prt_in  = new_prt_spec (prt_in), &
            prt_out = new_prt_spec (prt_out), &
            method  = var_str ("test_decay"), &
            variant = def)
       call lib%append (entry)
    end if

    call lib%configure (os_data)
    call lib%load (os_data)
  end subroutine prc_test_create_library


end module prc_test
