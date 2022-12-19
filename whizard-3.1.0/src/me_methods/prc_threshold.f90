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
module prc_threshold

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use physics_defs
  use diagnostics
  use os_interface
  use lorentz
  use interactions
  use model_data
  use variables, only: var_list_t

  use prclib_interfaces
  use process_libraries
  use prc_core_def
  use prc_core
  use prc_external

  implicit none
  private

  public :: threshold_writer_t
  public :: threshold_set_process_mode
  public :: threshold_get_amp_squared
  public :: threshold_olp_eval2
  public :: threshold_init
  public :: threshold_start_openloops
  public :: threshold_driver_t
  public :: threshold_def_t
  public :: threshold_state_t
  public :: prc_threshold_t

  interface
    subroutine threshold_set_process_mode (mode) bind(C)
      import
      integer(kind = c_int), intent(in) :: mode
    end subroutine threshold_set_process_mode
  end interface

  interface
    subroutine threshold_get_amp_squared (amp2, p_ofs, p_ons, leg, n_tot, sel_hel_beam) bind(C)
      import
      real(c_default_float), intent(out) :: amp2
      real(c_default_float), dimension(0:3,*), intent(in) :: p_ofs
      real(c_default_float), dimension(0:3,*), intent(in) :: p_ons
      integer(kind = c_int) :: n_tot, leg, sel_hel_beam
    end subroutine threshold_get_amp_squared
  end interface

  interface
    subroutine threshold_olp_eval2 (i_flv, alpha_s_c, parray, mu_c, &
           sel_hel_beam, sqme_c, acc_c) bind(C)
      import
      integer(c_int), intent(in) :: i_flv
      real(c_default_float), intent(in) :: alpha_s_c
      real(c_default_float), dimension(0:3,*), intent(in) :: parray
      real(c_default_float), intent(in) :: mu_c
      integer, intent(in) :: sel_hel_beam
      real(c_default_float), dimension(4), intent(out) :: sqme_c
      real(c_default_float), intent(out) :: acc_c
    end subroutine threshold_olp_eval2
  end interface

  interface
   subroutine threshold_init (par, scheme) bind(C)
      import
      real(c_default_float), dimension(*), intent(in) :: par
      integer(c_int), intent(in) :: scheme
    end subroutine threshold_init
  end interface

  interface
   subroutine threshold_start_openloops () bind(C)
      import
    end subroutine threshold_start_openloops
  end interface


  type, extends (prc_external_writer_t) :: threshold_writer_t
     integer :: nlo_type
  contains
    procedure :: init => threshold_writer_init
    procedure :: write_makefile_extra => threshold_writer_write_makefile_extra
    procedure :: write_makefile_code => threshold_writer_write_makefile_code
    procedure, nopass :: type_name => threshold_writer_type_name
  end type threshold_writer_t

  type, extends (prc_external_driver_t) :: threshold_driver_t
    procedure(threshold_olp_eval2), nopass, pointer :: &
         olp_eval2 => null ()
    procedure(threshold_set_process_mode), nopass, pointer :: &
         set_process_mode => null ()
    procedure(threshold_get_amp_squared), nopass, pointer :: &
         get_amp_squared => null ()
    procedure(threshold_start_openloops), nopass, pointer :: &
         start_openloops => null ()
    procedure(threshold_init), nopass, pointer :: &
         init => null ()
    type(string_t) :: id
    integer :: nlo_type = BORN
  contains
    procedure, nopass :: type_name => threshold_driver_type_name
    procedure :: load => threshold_driver_load
  end type threshold_driver_t

  type, extends (prc_external_def_t) :: threshold_def_t
     integer :: nlo_type
  contains
    procedure :: init => threshold_def_init
    procedure, nopass :: type_string => threshold_def_type_string
    procedure :: write => threshold_def_write
    procedure :: read => threshold_def_read
    procedure :: allocate_driver => threshold_def_allocate_driver
    procedure :: connect => threshold_def_connect
  end type threshold_def_t

  type, extends (prc_external_state_t) :: threshold_state_t
  contains
    procedure :: write => threshold_state_write
  end type threshold_state_t

  type, extends (prc_external_t) :: prc_threshold_t
     real(default), dimension(:,:), allocatable :: parray_ofs
     real(default), dimension(:,:), allocatable :: parray_ons
     integer :: leg
     logical :: has_beam_pol = .false.
  contains
    procedure :: write => prc_threshold_write
    procedure :: write_name => prc_threshold_write_name
    procedure :: set_beam_pol => prc_threshold_set_beam_pol
    procedure :: compute_amplitude => prc_threshold_compute_amplitude
    procedure :: allocate_workspace => prc_threshold_allocate_workspace
    procedure :: set_offshell_momenta => prc_threshold_set_offshell_momenta
    procedure :: set_onshell_momenta => prc_threshold_set_onshell_momenta
    procedure :: set_leg => prc_threshold_set_leg
    procedure :: set_process_mode => prc_threshold_set_process_mode
    procedure :: compute_sqme => prc_threshold_compute_sqme
    procedure :: compute_sqme_virt => prc_threshold_compute_sqme_virt
    procedure :: init => prc_threshold_init
    procedure :: activate_parameters => prc_threshold_activate_parameters
    procedure :: prepare_external_code => &
         prc_threshold_prepare_external_code
    procedure :: includes_polarization => prc_threshold_includes_polarization
  end type prc_threshold_t


  interface
    pure module subroutine threshold_writer_init &
         (writer, model_name, prt_in, prt_out, restrictions)
      class(threshold_writer_t), intent(inout) :: writer
      type(string_t), intent(in) :: model_name
      type(string_t), dimension(:), intent(in) :: prt_in, prt_out
      type(string_t), intent(in), optional :: restrictions
    end subroutine threshold_writer_init
    module subroutine threshold_writer_write_makefile_extra &
         (writer, unit, id, os_data, verbose, nlo_type)
      class(threshold_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: verbose
      integer, intent(in) :: nlo_type
    end subroutine threshold_writer_write_makefile_extra
    module subroutine threshold_writer_write_makefile_code &
         (writer, unit, id, os_data, verbose, testflag)
      class(threshold_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: verbose
      logical, intent(in), optional :: testflag
    end subroutine threshold_writer_write_makefile_code
    module function threshold_writer_type_name () result (string)
      type(string_t) :: string
    end function threshold_writer_type_name
    module function threshold_driver_type_name () result (type)
      type(string_t) :: type
    end function threshold_driver_type_name
    module subroutine threshold_driver_load (threshold_driver, dlaccess)
      class(threshold_driver_t), intent(inout) :: threshold_driver
      type(dlaccess_t), intent(inout) :: dlaccess
    end subroutine threshold_driver_load
    module function threshold_def_type_string () result (string)
      type(string_t) :: string
    end function threshold_def_type_string
    module subroutine threshold_def_write (object, unit)
      class(threshold_def_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine threshold_def_write
    module subroutine threshold_def_read (object, unit)
      class(threshold_def_t), intent(out) :: object
      integer, intent(in) :: unit
    end subroutine threshold_def_read
    module subroutine threshold_def_connect (def, lib_driver, i, proc_driver)
      class(threshold_def_t), intent(in) :: def
      class(prclib_driver_t), intent(in) :: lib_driver
      integer, intent(in) :: i
      class(prc_core_driver_t), intent(inout) :: proc_driver
    end subroutine threshold_def_connect
    module subroutine threshold_state_write (object, unit)
      class(threshold_state_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine threshold_state_write
    module subroutine prc_threshold_write (object, unit)
      class(prc_threshold_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_threshold_write
    module subroutine prc_threshold_write_name (object, unit)
      class(prc_threshold_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_threshold_write_name
    module subroutine prc_threshold_set_beam_pol (object, has_beam_pol)
      class(prc_threshold_t), intent(inout) :: object
      logical, intent(in), optional :: has_beam_pol
    end subroutine prc_threshold_set_beam_pol
    module function prc_threshold_compute_amplitude &
         (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
         core_state)  result (amp)
      class(prc_threshold_t), intent(in) :: object
      integer, intent(in) :: j
      type(vector4_t), dimension(:), intent(in) :: p
      integer, intent(in) :: f, h, c
      real(default), intent(in) :: fac_scale, ren_scale
      real(default), intent(in), allocatable :: alpha_qcd_forced
      class(prc_core_state_t), intent(inout), allocatable, optional :: &
           core_state
      complex(default) :: amp
    end function prc_threshold_compute_amplitude
    module subroutine prc_threshold_set_offshell_momenta (object, p)
      class(prc_threshold_t), intent(inout) :: object
      type(vector4_t), intent(in), dimension(:) :: p
    end subroutine prc_threshold_set_offshell_momenta
    module subroutine prc_threshold_set_onshell_momenta (object, p)
      class(prc_threshold_t), intent(inout) :: object
      type(vector4_t), intent(in), dimension(:) :: p
    end subroutine prc_threshold_set_onshell_momenta
    module subroutine prc_threshold_set_leg (object, leg)
      class(prc_threshold_t), intent(inout) :: object
      integer, intent(in) :: leg
    end subroutine prc_threshold_set_leg
    module subroutine prc_threshold_set_process_mode (object, mode)
      class(prc_threshold_t), intent(in) :: object
      integer(kind = c_int), intent(in) :: mode
    end subroutine prc_threshold_set_process_mode
    module subroutine prc_threshold_compute_sqme (object, i_flv, i_hel, p, &
           ren_scale, sqme, bad_point)
      class(prc_threshold_t), intent(in) :: object
      integer, intent(in) :: i_flv, i_hel
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(in) :: ren_scale
      real(default), intent(out) :: sqme
      logical, intent(out) :: bad_point
    end subroutine prc_threshold_compute_sqme
    module subroutine prc_threshold_compute_sqme_virt (object, i_flv, i_hel, &
           p, ren_scale, es_scale, loop_method, sqme, bad_point)
      class(prc_threshold_t), intent(in) :: object
      integer, intent(in) :: i_flv, i_hel
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: ren_scale, es_scale
      integer, intent(in) :: loop_method
      real(default), dimension(4), intent(out) :: sqme
      real(c_default_float), dimension(:,:), allocatable, save :: parray
      logical, intent(out) :: bad_point
    end subroutine prc_threshold_compute_sqme_virt
    module subroutine prc_threshold_init (object, def, lib, id, i_component)
      class(prc_threshold_t), intent(inout) :: object
      class(prc_core_def_t), intent(in), target :: def
      type(process_library_t), intent(in), target :: lib
      type(string_t), intent(in) :: id
      integer, intent(in) :: i_component
    end subroutine prc_threshold_init
    module subroutine prc_threshold_activate_parameters (object)
      class (prc_threshold_t), intent(inout) :: object
    end subroutine prc_threshold_activate_parameters
    module subroutine prc_threshold_prepare_external_code &
         (core, flv_states, var_list, os_data, libname, model, i_core, is_nlo)
      class(prc_threshold_t), intent(inout) :: core
      integer, intent(in), dimension(:,:), allocatable :: flv_states
      type(var_list_t), intent(in) :: var_list
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in) :: libname
      type(model_data_t), intent(in), target :: model
      integer, intent(in) :: i_core
      logical, intent(in) :: is_nlo
    end subroutine prc_threshold_prepare_external_code
    module function prc_threshold_includes_polarization &
         (object) result (polarized)
      logical :: polarized
      class(prc_threshold_t), intent(in) :: object
    end function prc_threshold_includes_polarization
  end interface

contains

  subroutine threshold_def_init (object, basename, model_name, &
       prt_in, prt_out, nlo_type, restrictions)
    class(threshold_def_t), intent(inout) :: object
    type(string_t), intent(in) :: basename, model_name
    type(string_t), dimension(:), intent(in) :: prt_in, prt_out
    integer, intent(in) :: nlo_type
    type(string_t), intent(in), optional :: restrictions
    if (debug_on) call msg_debug (D_ME_METHODS, "threshold_def_init")
    object%basename = basename
    object%nlo_type = nlo_type
    allocate (threshold_writer_t :: object%writer)
    select type (writer => object%writer)
    type is (threshold_writer_t)
       call writer%init (model_name, prt_in, prt_out, restrictions)
       writer%nlo_type = nlo_type
    end select
  end subroutine threshold_def_init

  subroutine threshold_def_allocate_driver (object, driver, basename)
    class(threshold_def_t), intent(in) :: object
    class(prc_core_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    if (debug_on) call msg_debug (D_ME_METHODS, "threshold_def_allocate_driver")
    if (.not. allocated (driver)) allocate (threshold_driver_t :: driver)
    select type (driver)
    type is (threshold_driver_t)
       driver%id = basename
       driver%nlo_type = object%nlo_type
    end select
  end subroutine threshold_def_allocate_driver

  subroutine prc_threshold_allocate_workspace (object, core_state)
    class(prc_threshold_t), intent(in) :: object
    class(prc_core_state_t), intent(inout), allocatable :: core_state
    allocate (threshold_state_t :: core_state)
  end subroutine prc_threshold_allocate_workspace


end module prc_threshold
