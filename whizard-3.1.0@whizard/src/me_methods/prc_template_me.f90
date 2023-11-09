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

module prc_template_me

  use, intrinsic ::  iso_c_binding !NODEP!

  use kinds
  use iso_varying_string, string_t => varying_string
  use os_interface
  use lorentz
  use interactions
  use model_data

  use particle_specifiers, only: new_prt_spec
  use process_constants
  use prclib_interfaces
  use prc_core_def
  use process_libraries
  use prc_core

  implicit none
  private

  public :: template_me_def_t
  public :: template_me_driver_t
  public :: template_me_make_process_component
  public :: prc_template_me_t

  type, extends (prc_core_def_t) :: template_me_def_t
   contains
     procedure, nopass :: type_string => template_me_def_type_string
     procedure :: init => template_me_def_init
     procedure :: write => template_me_def_write
     procedure :: read => template_me_def_read
     procedure :: allocate_driver => template_me_def_allocate_driver
     procedure, nopass :: needs_code => template_me_def_needs_code
     procedure, nopass :: get_features => template_me_def_get_features
     procedure :: connect => template_me_def_connect
  end type template_me_def_t

  type, extends (prc_writer_f_module_t) :: template_me_writer_t
     class(model_data_t), pointer :: model => null ()
     type(string_t) :: model_name
     logical :: unity
     type(string_t), dimension(:), allocatable :: prt_in
     type(string_t), dimension(:), allocatable :: prt_out
     integer :: n_in
     integer :: n_out
     integer :: n_tot
   contains
     procedure, nopass :: type_name => template_me_writer_type_name
     procedure, nopass :: get_module_name => template_me_writer_get_module_name
     procedure :: write => template_me_writer_write
     procedure :: init => template_me_writer_init
     procedure :: write_makefile_code => template_me_write_makefile_code
     procedure :: write_source_code => template_me_write_source_code
     procedure :: before_compile => template_me_before_compile
     procedure :: after_compile => template_me_after_compile
     procedure, nopass :: get_procname => template_me_writer_get_procname
     procedure :: write_interface => template_me_write_interface
     procedure :: write_wrapper => template_me_write_wrapper
  end type template_me_writer_t

  type, extends (prc_core_driver_t) :: template_me_driver_t
     procedure(init_t), nopass, pointer :: &
          init => null ()
     procedure(update_alpha_s_t), nopass, pointer :: &
          update_alpha_s => null ()
     procedure(is_allowed_t), nopass, pointer :: &
          is_allowed => null ()
     procedure(new_event_t), nopass, pointer :: &
          new_event => null ()
     procedure(get_amplitude_t), nopass, pointer :: &
          get_amplitude => null ()
   contains
     procedure, nopass :: type_name => template_me_driver_type_name
  end type template_me_driver_t

  type, extends (prc_core_t) :: prc_template_me_t
     real(default), dimension(:), allocatable :: par
     integer :: scheme = 0
   contains
     procedure :: allocate_workspace => prc_template_me_allocate_workspace
     procedure :: write => prc_template_me_write
     procedure :: write_name => prc_template_me_write_name
     procedure :: set_parameters => prc_template_me_set_parameters
     procedure :: init => prc_template_me_init
     procedure :: activate_parameters => prc_template_me_activate_parameters
     procedure :: is_allowed => prc_template_me_is_allowed
     procedure :: compute_hard_kinematics => &
          prc_template_me_compute_hard_kinematics
     procedure :: compute_eff_kinematics => &
          prc_template_me_compute_eff_kinematics
     procedure :: compute_amplitude => prc_template_me_compute_amplitude
  end type prc_template_me_t

  type, extends (prc_core_state_t) :: template_me_state_t
     logical :: new_kinematics = .true.
     real(default) :: alpha_qcd = -1
     real(default) :: alpha_qed = -1
   contains
     procedure :: write => template_me_state_write
     procedure :: reset_new_kinematics => template_me_state_reset_new_kinematics
  end type template_me_state_t


  abstract interface
     subroutine init_t (par, scheme) bind(C)
       import
       real(c_default_float), dimension(*), intent(in) :: par
       integer(c_int), intent(in) :: scheme
     end subroutine init_t
  end interface

  abstract interface
     subroutine update_alpha_s_t (alpha_s) bind(C)
       import
       real(c_default_float), intent(in) :: alpha_s
     end subroutine update_alpha_s_t
  end interface

  abstract interface
     subroutine is_allowed_t (flv, hel, col, flag) bind(C)
       import
       integer(c_int), intent(in) :: flv, hel, col
       logical(c_bool), intent(out) :: flag
     end subroutine is_allowed_t
  end interface

  abstract interface
     subroutine new_event_t (p) bind(C)
       import
       real(c_default_float), dimension(0:3,*), intent(in) :: p
     end subroutine new_event_t
  end interface

  abstract interface
     subroutine get_amplitude_t (flv, hel, col, amp) bind(C)
       import
       integer(c_int), intent(in) :: flv, hel, col
       complex(c_default_complex), intent(out):: amp
     end subroutine get_amplitude_t
  end interface


  interface
    module function template_me_def_type_string () result (string)
      type(string_t) :: string
    end function template_me_def_type_string
    module subroutine template_me_def_write (object, unit)
      class(template_me_def_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine template_me_def_write
    module subroutine template_me_def_read (object, unit)
      class(template_me_def_t), intent(out) :: object
      integer, intent(in) :: unit
    end subroutine template_me_def_read
    module function template_me_def_needs_code () result (flag)
      logical :: flag
    end function template_me_def_needs_code
    module subroutine template_me_def_get_features (features)
      type(string_t), dimension(:), allocatable, intent(out) :: features
    end subroutine template_me_def_get_features
    module subroutine template_me_def_connect (def, lib_driver, i, proc_driver)
      class(template_me_def_t), intent(in) :: def
      class(prclib_driver_t), intent(in) :: lib_driver
      integer, intent(in) :: i
      class(prc_core_driver_t), intent(inout) :: proc_driver
    end subroutine template_me_def_connect
    module function template_me_writer_type_name () result (string)
      type(string_t) :: string
    end function template_me_writer_type_name
    module function template_me_writer_get_module_name (id) result (name)
      type(string_t) :: name
      type(string_t), intent(in) :: id
    end function template_me_writer_get_module_name
    module subroutine template_me_writer_write (object, unit)
      class(template_me_writer_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine template_me_writer_write
    module subroutine template_me_writer_init (writer, model, &
         prt_in, prt_out, unity)
      class(template_me_writer_t), intent(out) :: writer
      class(model_data_t), intent(in), target :: model
      type(string_t), dimension(:), intent(in) :: prt_in
      type(string_t), dimension(:), intent(in) :: prt_out
      logical, intent(in) :: unity
    end subroutine template_me_writer_init
    module subroutine template_me_write_makefile_code &
         (writer, unit, id, os_data, verbose, testflag)
      class(template_me_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: verbose
      logical, intent(in), optional :: testflag
    end subroutine template_me_write_makefile_code
    module subroutine template_me_write_source_code (writer, id)
      class(template_me_writer_t), intent(in) :: writer
      type(string_t), intent(in) :: id
    end subroutine template_me_write_source_code
    module subroutine template_me_before_compile (writer, id)
      class(template_me_writer_t), intent(in) :: writer
      type(string_t), intent(in) :: id
    end subroutine template_me_before_compile
    module subroutine template_me_after_compile (writer, id)
      class(template_me_writer_t), intent(in) :: writer
      type(string_t), intent(in) :: id
    end subroutine template_me_after_compile
    module function template_me_writer_get_procname (feature) result (name)
      type(string_t) :: name
      type(string_t), intent(in) :: feature
    end function template_me_writer_get_procname
    module subroutine template_me_write_interface (writer, unit, id, feature)
      class(template_me_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
      type(string_t), intent(in) :: feature
    end subroutine template_me_write_interface
    module subroutine template_me_write_wrapper (writer, unit, id, feature)
      class(template_me_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id, feature
    end subroutine template_me_write_wrapper
    module function template_me_driver_type_name () result (string)
      type(string_t) :: string
    end function template_me_driver_type_name
    module subroutine template_me_state_write (object, unit)
      class(template_me_state_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine template_me_state_write  
    module subroutine template_me_state_reset_new_kinematics (object)
      class(template_me_state_t), intent(inout) :: object
    end subroutine template_me_state_reset_new_kinematics
    module subroutine prc_template_me_write (object, unit)
      class(prc_template_me_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_template_me_write
    module subroutine prc_template_me_write_name (object, unit)
      class(prc_template_me_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_template_me_write_name
    module subroutine prc_template_me_set_parameters (prc_template_me, model)
      class(prc_template_me_t), intent(inout) :: prc_template_me
      class(model_data_t), intent(in), target, optional :: model
    end subroutine prc_template_me_set_parameters
    module subroutine prc_template_me_init (object, def, lib, id, i_component)
      class(prc_template_me_t), intent(inout) :: object
      class(prc_core_def_t), intent(in), target :: def
      type(process_library_t), intent(in), target :: lib
      type(string_t), intent(in) :: id
      integer, intent(in) :: i_component
    end subroutine prc_template_me_init
    module subroutine prc_template_me_activate_parameters (object)
      class (prc_template_me_t), intent(inout) :: object
    end subroutine prc_template_me_activate_parameters
    module function prc_template_me_is_allowed &
        (object, i_term, f, h, c) result (flag)
      class(prc_template_me_t), intent(in) :: object
      integer, intent(in) :: i_term, f, h, c
      logical :: flag
    end function prc_template_me_is_allowed
    module subroutine prc_template_me_compute_hard_kinematics &
         (object, p_seed, i_term, int_hard, core_state)
      class(prc_template_me_t), intent(in) :: object
      type(vector4_t), dimension(:), intent(in) :: p_seed
      integer, intent(in) :: i_term
      type(interaction_t), intent(inout) :: int_hard
      class(prc_core_state_t), intent(inout), allocatable :: core_state
    end subroutine prc_template_me_compute_hard_kinematics
    module subroutine prc_template_me_compute_eff_kinematics &
         (object, i_term, int_hard, int_eff, core_state)
      class(prc_template_me_t), intent(in) :: object
      integer, intent(in) :: i_term
      type(interaction_t), intent(in) :: int_hard
      type(interaction_t), intent(inout) :: int_eff
      class(prc_core_state_t), intent(inout), allocatable :: core_state
    end subroutine prc_template_me_compute_eff_kinematics
    module function prc_template_me_compute_amplitude &
         (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
         core_state)  result (amp)
      class(prc_template_me_t), intent(in) :: object
      integer, intent(in) :: j
      type(vector4_t), dimension(:), intent(in) :: p
      integer, intent(in) :: f, h, c
      real(default), intent(in) :: fac_scale, ren_scale
      real(default), intent(in), allocatable :: alpha_qcd_forced
      class(prc_core_state_t), intent(inout), allocatable, optional :: &
           core_state
      complex(default) :: amp
    end function prc_template_me_compute_amplitude
  end interface

contains

  subroutine template_me_def_init &
       (object, model, prt_in, prt_out, unity)
    class(template_me_def_t), intent(out) :: object
    class(model_data_t), intent(in), target :: model
    type(string_t), dimension(:), intent(in) :: prt_in
    type(string_t), dimension(:), intent(in) :: prt_out
    logical, intent(in) :: unity
    allocate (template_me_writer_t :: object%writer)
    select type (writer => object%writer)
    type is (template_me_writer_t)
       call writer%init (model, prt_in, prt_out, unity)
    end select
  end subroutine template_me_def_init

  subroutine template_me_def_allocate_driver (object, driver, basename)
    class(template_me_def_t), intent(in) :: object
    class(prc_core_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    allocate (template_me_driver_t :: driver)
  end subroutine template_me_def_allocate_driver

  subroutine template_me_make_process_component (entry, component_index, &
         model, model_name, prt_in, prt_out, unity)
    class(process_def_entry_t), intent(inout) :: entry
    integer, intent(in) :: component_index
    type(string_t), intent(in) :: model_name
    class(model_data_t), intent(in), target :: model
    type(string_t), dimension(:), intent(in) :: prt_in
    type(string_t), dimension(:), intent(in) :: prt_out
    logical, intent(in) :: unity
    class(prc_core_def_t), allocatable :: def
    allocate (template_me_def_t :: def)
    select type (def)
    type is (template_me_def_t)
       call def%init (model, prt_in, prt_out, unity)
    end select
    call entry%import_component (component_index, &
         n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method = var_str ("template"), &
         variant = def)
  end subroutine template_me_make_process_component

  subroutine prc_template_me_allocate_workspace (object, core_state)
    class(prc_template_me_t), intent(in) :: object
    class(prc_core_state_t), intent(inout), allocatable :: core_state
    allocate (template_me_state_t :: core_state)
  end subroutine prc_template_me_allocate_workspace


end module prc_template_me
