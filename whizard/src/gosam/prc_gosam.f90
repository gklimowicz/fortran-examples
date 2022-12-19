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

module prc_gosam

  use, intrinsic :: iso_c_binding !NODEP!
  use, intrinsic :: iso_fortran_env

  use kinds
  use iso_varying_string, string_t => varying_string
  use physics_defs
  use os_interface
  use lorentz
  use interactions
  use model_data
  use variables
  use prc_core_def
  use prc_core
  use blha_config
  use blha_olp_interfaces

  implicit none
  private



  public :: gosam_def_t
  public :: prc_gosam_t

  type, extends (prc_blha_writer_t) :: gosam_writer_t
    type(string_t) :: gosam_dir
    type(string_t) :: golem_dir
    type(string_t) :: samurai_dir
    type(string_t) :: ninja_dir
    type(string_t) :: form_dir
    type(string_t) :: qgraf_dir
    type(string_t) :: filter_lo, filter_nlo
    type(string_t) :: symmetries
    integer :: form_threads
    integer :: form_workspace
    type(string_t) :: fc
  contains
    procedure :: write_config => gosam_writer_write_config
    procedure, nopass :: type_name => gosam_writer_type_name
    procedure :: init => gosam_writer_init
    procedure :: generate_configuration_file => &
              gosam_writer_generate_configuration_file
  end type gosam_writer_t

  type, extends (blha_def_t) :: gosam_def_t
    logical :: execute_olp = .true.
  contains
    procedure :: init => gosam_def_init
    procedure, nopass :: type_string => gosam_def_type_string
    procedure :: write => gosam_def_write
    procedure :: read => gosam_def_read
    procedure :: allocate_driver => gosam_def_allocate_driver
  end type gosam_def_t

  type, extends (blha_driver_t) :: gosam_driver_t
    type(string_t) :: gosam_dir
    type(string_t) :: olp_file
    type(string_t) :: olc_file
    type(string_t) :: olp_dir
    type(string_t) :: olp_lib
  contains
    procedure, nopass :: type_name => gosam_driver_type_name
    procedure :: init_gosam => gosam_driver_init_gosam
    procedure :: init_dlaccess_to_library => gosam_driver_init_dlaccess_to_library
    procedure :: write_makefile => gosam_driver_write_makefile
    procedure :: set_alpha_s => gosam_driver_set_alpha_s
    procedure :: set_alpha_qed => gosam_driver_set_alpha_qed
    procedure :: set_GF => gosam_driver_set_GF
    procedure :: set_weinberg_angle => gosam_driver_set_weinberg_angle
    procedure :: print_alpha_s => gosam_driver_print_alpha_s
  end type gosam_driver_t

  type, extends (prc_blha_t) :: prc_gosam_t
    logical :: initialized = .false.
  contains
    procedure :: prepare_library => prc_gosam_prepare_library
    procedure :: prepare_external_code => &
         prc_gosam_prepare_external_code
    procedure :: write_makefile => prc_gosam_write_makefile
    procedure :: execute_makefile => prc_gosam_execute_makefile
    procedure :: create_olp_library => prc_gosam_create_olp_library
    procedure :: load_driver => prc_gosam_load_driver
    procedure :: start => prc_gosam_start
    procedure :: write => prc_gosam_write
    procedure :: write_name => prc_gosam_write_name
    procedure :: init_driver => prc_gosam_init_driver
    procedure :: set_initialized => prc_gosam_set_initialized
    procedure :: compute_sqme_spin_c => prc_gosam_compute_sqme_spin_c
    procedure :: allocate_workspace => prc_gosam_allocate_workspace
    procedure :: set_particle_properties => prc_gosam_set_particle_properties
  end type prc_gosam_t

  type, extends (blha_state_t) :: gosam_state_t
  contains
    procedure :: write => gosam_state_write
  end type gosam_state_t




  interface
    module subroutine gosam_writer_write_config (gosam_writer)
      class(gosam_writer_t), intent(in) :: gosam_writer
    end subroutine gosam_writer_write_config
    module function gosam_def_type_string () result (string)
      type(string_t) :: string
    end function gosam_def_type_string
    module subroutine gosam_def_write (object, unit)
      class(gosam_def_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine gosam_def_write
    module subroutine gosam_def_read (object, unit)
      class(gosam_def_t), intent(out) :: object
      integer, intent(in) :: unit
    end subroutine gosam_def_read
    module function gosam_writer_type_name () result (string)
      type(string_t) :: string
    end function gosam_writer_type_name
    pure module subroutine gosam_writer_init &
         (writer, model_name, prt_in, prt_out, restrictions)
      class(gosam_writer_t), intent(inout) :: writer
      type(string_t), intent(in) :: model_name
      type(string_t), dimension(:), intent(in) :: prt_in, prt_out
      type(string_t), intent(in), optional :: restrictions
    end subroutine gosam_writer_init
    module function gosam_driver_type_name () result (string)
      type(string_t) :: string
    end function gosam_driver_type_name
    module subroutine gosam_driver_init_gosam (object, os_data, olp_file, &
         olc_file, olp_dir, olp_lib)
      class(gosam_driver_t), intent(inout) :: object
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in) :: olp_file, olc_file, olp_dir, olp_lib
    end subroutine gosam_driver_init_gosam
    module subroutine gosam_driver_init_dlaccess_to_library &
       (object, os_data, dlaccess, success)
      class(gosam_driver_t), intent(in) :: object
      type(os_data_t), intent(in) :: os_data
      type(dlaccess_t), intent(out) :: dlaccess
      logical, intent(out) :: success
    end subroutine gosam_driver_init_dlaccess_to_library
    module subroutine gosam_writer_generate_configuration_file &
         (object, unit)
      class(gosam_writer_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine gosam_writer_generate_configuration_file
    module subroutine gosam_driver_write_makefile (object, unit, libname)
      class(gosam_driver_t), intent(in) :: object
      integer, intent(in) :: unit
      type(string_t), intent(in) :: libname
    end subroutine gosam_driver_write_makefile
    module subroutine gosam_driver_set_alpha_s (driver, alpha_s)
       class(gosam_driver_t), intent(in) :: driver
       real(default), intent(in) :: alpha_s
    end subroutine gosam_driver_set_alpha_s
    module subroutine gosam_driver_set_alpha_qed (driver, alpha)
      class(gosam_driver_t), intent(inout) :: driver
      real(default), intent(in) :: alpha
    end subroutine gosam_driver_set_alpha_qed
    module subroutine gosam_driver_set_GF (driver, GF)
      class(gosam_driver_t), intent(inout) :: driver
      real(default), intent(in) :: GF
    end subroutine gosam_driver_set_GF
    module subroutine gosam_driver_set_weinberg_angle (driver, sw2)
      class(gosam_driver_t), intent(inout) :: driver
      real(default), intent(in) :: sw2
    end subroutine gosam_driver_set_weinberg_angle
    module subroutine gosam_driver_print_alpha_s (object)
      class(gosam_driver_t), intent(in) :: object
    end subroutine gosam_driver_print_alpha_s
    module subroutine prc_gosam_prepare_library (object, os_data, libname)
      class(prc_gosam_t), intent(inout) :: object
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in) :: libname
    end subroutine prc_gosam_prepare_library
    module subroutine prc_gosam_prepare_external_code &
         (core, flv_states, var_list, os_data, libname, model, i_core, is_nlo)
      class(prc_gosam_t), intent(inout) :: core
      integer, intent(in), dimension(:,:), allocatable :: flv_states
      type(var_list_t), intent(in) :: var_list
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in) :: libname
      type(model_data_t), intent(in), target :: model
      integer, intent(in) :: i_core
      logical, intent(in) :: is_nlo
    end subroutine prc_gosam_prepare_external_code
    module subroutine prc_gosam_write_makefile (object, unit, libname)
      class(prc_gosam_t), intent(in) :: object
      integer, intent(in) :: unit
      type(string_t), intent(in) :: libname
    end subroutine prc_gosam_write_makefile
    module subroutine prc_gosam_execute_makefile (object, libname)
      class(prc_gosam_t), intent(in) :: object
      type(string_t), intent(in) :: libname
    end subroutine prc_gosam_execute_makefile
    module subroutine prc_gosam_create_olp_library (object, libname)
      class(prc_gosam_t), intent(inout) :: object
      type(string_t), intent(in) :: libname
    end subroutine prc_gosam_create_olp_library
    module subroutine prc_gosam_load_driver (object, os_data)
      class(prc_gosam_t), intent(inout) :: object
      type(os_data_t), intent(in) :: os_data
    end subroutine prc_gosam_load_driver
    module subroutine prc_gosam_start (object)
      class(prc_gosam_t), intent(inout) :: object
    end subroutine prc_gosam_start
    module subroutine prc_gosam_write (object, unit)
      class(prc_gosam_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_gosam_write
    module subroutine prc_gosam_write_name (object, unit)
      class(prc_gosam_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_gosam_write_name
    module subroutine prc_gosam_init_driver (object, os_data)
      class(prc_gosam_t), intent(inout) :: object
      type(os_data_t), intent(in) :: os_data
    end subroutine prc_gosam_init_driver
    module subroutine prc_gosam_set_initialized (prc_gosam)
      class(prc_gosam_t), intent(inout) :: prc_gosam
    end subroutine prc_gosam_set_initialized
    module subroutine prc_gosam_compute_sqme_spin_c (object, &
         i_flv, i_hel, em, p, ren_scale, me_sc, bad_point)
      class(prc_gosam_t), intent(inout) :: object
      integer, intent(in) :: i_flv, i_hel
      integer, intent(in) :: em
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(in) :: ren_scale
      complex(default), intent(out) :: me_sc
      logical, intent(out) :: bad_point
    end subroutine prc_gosam_compute_sqme_spin_c
    module subroutine gosam_state_write (object, unit)
      class(gosam_state_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine gosam_state_write
    module subroutine prc_gosam_set_particle_properties (object, model)
      class(prc_gosam_t), intent(inout) :: object
      class(model_data_t), intent(in), target :: model
    end subroutine prc_gosam_set_particle_properties
  end interface

contains

  subroutine gosam_def_init (object, basename, model_name, &
     prt_in, prt_out, nlo_type, restrictions, var_list)
    class(gosam_def_t), intent(inout) :: object
    type(string_t), intent(in) :: basename
    type(string_t), intent(in) :: model_name
    type(string_t), dimension(:), intent(in) :: prt_in, prt_out
    integer, intent(in) :: nlo_type
    type(string_t), intent(in), optional :: restrictions
    type(var_list_t), intent(in) :: var_list
    object%basename = basename
    allocate (gosam_writer_t :: object%writer)
    select case (nlo_type)
    case (BORN)
       object%suffix = '_BORN'
    case (NLO_REAL)
       object%suffix = '_REAL'
    case (NLO_VIRTUAL)
       object%suffix = '_LOOP'
    case (NLO_SUBTRACTION)
       object%suffix = '_SUB'
    end select
    select type (writer => object%writer)
    type is (gosam_writer_t)
      call writer%init (model_name, prt_in, prt_out, restrictions)
      writer%filter_lo = var_list%get_sval (var_str ("$gosam_filter_lo"))
      writer%filter_nlo = var_list%get_sval (var_str ("$gosam_filter_nlo"))
      writer%symmetries = &
           var_list%get_sval (var_str ("$gosam_symmetries"))
      writer%form_threads = &
           var_list%get_ival (var_str ("form_threads"))
      writer%form_workspace = &
           var_list%get_ival (var_str ("form_workspace"))
      writer%fc = &
           var_list%get_sval (var_str ("$gosam_fc"))
    end select
  end subroutine gosam_def_init

  subroutine gosam_def_allocate_driver (object, driver, basename)
    class(gosam_def_t), intent(in) :: object
    class(prc_core_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    if (.not. allocated (driver)) allocate (gosam_driver_t :: driver)
  end subroutine gosam_def_allocate_driver

  subroutine prc_gosam_allocate_workspace (object, core_state)
    class(prc_gosam_t), intent(in) :: object
    class(prc_core_state_t), intent(inout), allocatable :: core_state
    allocate (gosam_state_t :: core_state)
  end subroutine prc_gosam_allocate_workspace


end module prc_gosam

