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

module prc_omega

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds
  use iso_varying_string, string_t => varying_string
  use constants, only: one
  use os_interface
  use diagnostics
  use lorentz
  use sm_qcd
  use sm_qed
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

  public :: omega_def_t
  public :: omega_driver_t
  public :: omega_make_process_component
  public :: prc_omega_t
  public :: omega_state_t

  type, extends (prc_core_def_t) :: omega_def_t
     logical :: ufo = .false.
     logical :: ovm = .false.
   contains
     procedure, nopass :: type_string => omega_def_type_string
     procedure :: init => omega_def_init
     procedure :: write => omega_def_write
     procedure :: read => omega_def_read
     procedure :: allocate_driver => omega_def_allocate_driver
     procedure, nopass :: needs_code => omega_def_needs_code
     procedure, nopass :: get_features => omega_def_get_features
     procedure :: connect => omega_def_connect
  end type omega_def_t

  type, extends (prc_writer_f_module_t), abstract :: omega_writer_t
     type(string_t) :: model_name
     type(string_t) :: process_mode
     type(string_t) :: process_string
     type(string_t) :: restrictions
     logical :: openmp_support = .false.
     logical :: report_progress = .false.
     logical :: diags = .false.
     logical :: diags_color = .false.
     logical :: complex_mass_scheme = .false.
     logical :: write_phs_output = .false.
     type(string_t) :: extra_options
   contains
     procedure, nopass :: get_module_name => omega_writer_get_module_name
     procedure :: write => omega_writer_write
     procedure :: init => omega_writer_init
     procedure :: write_makefile_code => omega_write_makefile_code
     procedure :: write_source_code => omega_write_source_code
     procedure :: before_compile => omega_before_compile
     procedure :: after_compile => omega_after_compile
     procedure, nopass :: get_procname => omega_writer_get_procname
     procedure :: write_interface => omega_write_interface
     procedure :: write_wrapper => omega_write_wrapper
  end type omega_writer_t

  type, extends (omega_writer_t) :: omega_omega_writer_t
   contains
     procedure, nopass :: type_name => omega_omega_writer_type_name
  end type omega_omega_writer_t

  type, extends (omega_omega_writer_t) :: omega_ufo_writer_t
     type(string_t) :: ufo_path
   contains
     procedure, nopass :: type_name => omega_ufo_writer_type_name
  end type omega_ufo_writer_t

  type, extends (omega_writer_t) :: omega_ovm_writer_t
   contains
     procedure, nopass :: type_name => omega_ovm_writer_type_name
  end type omega_ovm_writer_t

  type, extends (prc_core_driver_t) :: omega_driver_t
     procedure(init_t), nopass, pointer :: &
          init => null ()
     procedure(update_alpha_s_t), nopass, pointer :: &
          update_alpha_s => null ()
     procedure(reset_helicity_selection_t), nopass, pointer :: &
          reset_helicity_selection => null ()
     procedure(is_allowed_t), nopass, pointer :: &
          is_allowed => null ()
     procedure(new_event_t), nopass, pointer :: &
          new_event => null ()
     procedure(get_amplitude_t), nopass, pointer :: &
          get_amplitude => null ()
   contains
     procedure, nopass :: type_name => omega_driver_type_name
  end type omega_driver_t

  type, extends (prc_core_t) :: prc_omega_t
     real(default), dimension(:), allocatable :: par
     integer :: scheme = 0
     type(helicity_selection_t) :: helicity_selection
     type(qcd_t) :: qcd
     type(qed_t) :: qed
   contains
     procedure :: allocate_workspace => prc_omega_allocate_workspace
     procedure :: write => prc_omega_write
     procedure :: write_name => prc_omega_write_name
     procedure :: set_parameters => prc_omega_set_parameters
     procedure :: init => prc_omega_init
     procedure :: activate_parameters => prc_omega_activate_parameters
     procedure :: is_allowed => prc_omega_is_allowed
     procedure :: compute_hard_kinematics => prc_omega_compute_hard_kinematics
     procedure :: compute_eff_kinematics => prc_omega_compute_eff_kinematics
     procedure :: reset_helicity_selection => prc_omega_reset_helicity_selection
     procedure :: compute_amplitude => prc_omega_compute_amplitude
     procedure :: get_alpha_s => prc_omega_get_alpha_s
     procedure :: get_alpha_qed => prc_omega_get_alpha_qed
  end type prc_omega_t

  type, extends (prc_core_state_t) :: omega_state_t
     logical :: new_kinematics = .true.
     real(default) :: alpha_qcd = -1
     real(default) :: alpha_qed = -1
   contains
    procedure :: write => omega_state_write
    procedure :: reset_new_kinematics => omega_state_reset_new_kinematics
  end type omega_state_t


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
     subroutine reset_helicity_selection_t (threshold, cutoff) bind(C)
       import
       real(c_default_float), intent(in) :: threshold
       integer(c_int), intent(in) :: cutoff
     end subroutine reset_helicity_selection_t
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
    module function omega_def_type_string () result (string)
      type(string_t) :: string
    end function omega_def_type_string
    module subroutine omega_def_write (object, unit)
      class(omega_def_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine omega_def_write
    module subroutine omega_def_read (object, unit)
      class(omega_def_t), intent(out) :: object
      integer, intent(in) :: unit
    end subroutine omega_def_read
    module function omega_def_needs_code () result (flag)
      logical :: flag
    end function omega_def_needs_code
    module subroutine omega_def_get_features (features)
      type(string_t), dimension(:), allocatable, intent(out) :: features
    end subroutine omega_def_get_features
    module subroutine omega_def_connect (def, lib_driver, i, proc_driver)
      class(omega_def_t), intent(in) :: def
      class(prclib_driver_t), intent(in) :: lib_driver
      integer, intent(in) :: i
      class(prc_core_driver_t), intent(inout) :: proc_driver
    end subroutine omega_def_connect
    module function omega_omega_writer_type_name () result (string)
      type(string_t) :: string
    end function omega_omega_writer_type_name
    module function omega_ufo_writer_type_name () result (string)
      type(string_t) :: string
    end function omega_ufo_writer_type_name
    module function omega_ovm_writer_type_name () result (string)
      type(string_t) :: string
    end function omega_ovm_writer_type_name
    module function omega_writer_get_module_name (id) result (name)
      type(string_t) :: name
      type(string_t), intent(in) :: id
    end function omega_writer_get_module_name
    module subroutine omega_writer_write (object, unit)
      class(omega_writer_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine omega_writer_write
    module subroutine omega_writer_init (writer, model_name, prt_in, prt_out, &
         ufo_path, restrictions, cms_scheme, openmp_support, &
         report_progress, write_phs_output, extra_options, diags, diags_color)
      class(omega_writer_t), intent(out) :: writer
      type(string_t), intent(in) :: model_name
      type(string_t), dimension(:), intent(in) :: prt_in
      type(string_t), dimension(:), intent(in) :: prt_out
      type(string_t), intent(in), optional :: ufo_path
      type(string_t), intent(in), optional :: restrictions
      logical, intent(in), optional :: cms_scheme
      logical, intent(in), optional :: openmp_support
      logical, intent(in), optional :: report_progress
      logical, intent(in), optional :: write_phs_output
      type(string_t), intent(in), optional :: extra_options
      logical, intent(in), optional :: diags, diags_color
    end subroutine omega_writer_init
    module subroutine omega_write_makefile_code &
         (writer, unit, id, os_data, verbose, testflag)
      class(omega_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: verbose
      logical, intent(in), optional :: testflag
    end subroutine omega_write_makefile_code
    module subroutine omega_write_source_code (writer, id)
      class(omega_writer_t), intent(in) :: writer
      type(string_t), intent(in) :: id
    end subroutine omega_write_source_code
    module subroutine omega_before_compile (writer, id)
      class(omega_writer_t), intent(in) :: writer
      type(string_t), intent(in) :: id
    end subroutine omega_before_compile
    module subroutine omega_after_compile (writer, id)
      class(omega_writer_t), intent(in) :: writer
      type(string_t), intent(in) :: id
    end subroutine omega_after_compile
    module function omega_writer_get_procname (feature) result (name)
      type(string_t) :: name
      type(string_t), intent(in) :: feature
    end function omega_writer_get_procname
    module subroutine omega_write_interface (writer, unit, id, feature)
      class(omega_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
      type(string_t), intent(in) :: feature
    end subroutine omega_write_interface
    module subroutine omega_write_wrapper (writer, unit, id, feature)
      class(omega_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id, feature
    end subroutine omega_write_wrapper
    module function omega_driver_type_name () result (string)
      type(string_t) :: string
    end function omega_driver_type_name
    module subroutine omega_state_write (object, unit)
      class(omega_state_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine omega_state_write
    module subroutine omega_state_reset_new_kinematics (object)
      class(omega_state_t), intent(inout) :: object
    end subroutine omega_state_reset_new_kinematics
    module subroutine prc_omega_write (object, unit)
      class(prc_omega_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_omega_write
    module subroutine prc_omega_write_name (object, unit)
      class(prc_omega_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_omega_write_name
    module subroutine prc_omega_init (object, def, lib, id, i_component)
      class(prc_omega_t), intent(inout) :: object
      class(prc_core_def_t), intent(in), target :: def
      type(process_library_t), intent(in), target :: lib
      type(string_t), intent(in) :: id
      integer, intent(in) :: i_component
    end subroutine prc_omega_init
    module subroutine prc_omega_activate_parameters (object)
      class (prc_omega_t), intent(inout) :: object
    end subroutine prc_omega_activate_parameters
    module function prc_omega_is_allowed (object, i_term, f, h, c) result (flag)
      class(prc_omega_t), intent(in) :: object
      integer, intent(in) :: i_term, f, h, c
      logical :: flag
    end function prc_omega_is_allowed
    module subroutine prc_omega_compute_hard_kinematics &
         (object, p_seed, i_term, int_hard, core_state)
      class(prc_omega_t), intent(in) :: object
      type(vector4_t), dimension(:), intent(in) :: p_seed
      integer, intent(in) :: i_term
      type(interaction_t), intent(inout) :: int_hard
      class(prc_core_state_t), intent(inout), allocatable :: core_state
    end subroutine prc_omega_compute_hard_kinematics
    module subroutine prc_omega_compute_eff_kinematics &
         (object, i_term, int_hard, int_eff, core_state)
      class(prc_omega_t), intent(in) :: object
      integer, intent(in) :: i_term
      type(interaction_t), intent(in) :: int_hard
      type(interaction_t), intent(inout) :: int_eff
      class(prc_core_state_t), intent(inout), allocatable :: core_state
    end subroutine prc_omega_compute_eff_kinematics
    module subroutine prc_omega_reset_helicity_selection (object)
      class(prc_omega_t), intent(inout) :: object
    end subroutine prc_omega_reset_helicity_selection
    module function prc_omega_compute_amplitude &
         (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
         core_state)  result (amp)
      class(prc_omega_t), intent(in) :: object
      integer, intent(in) :: j
      type(vector4_t), dimension(:), intent(in) :: p
      integer, intent(in) :: f, h, c
      real(default), intent(in) :: fac_scale, ren_scale
      real(default), intent(in), allocatable :: alpha_qcd_forced
      class(prc_core_state_t), intent(inout), allocatable, optional :: &
           core_state
      complex(default) :: amp
    end function prc_omega_compute_amplitude
    module function prc_omega_get_alpha_s &
         (object, core_state) result (alpha_qcd)
      class(prc_omega_t), intent(in) :: object
      class(prc_core_state_t), intent(in), allocatable :: core_state
      real(default) :: alpha_qcd
    end function prc_omega_get_alpha_s
    module function prc_omega_get_alpha_qed &
         (object, core_state) result (alpha_qed)
      class(prc_omega_t), intent(in) :: object
      class(prc_core_state_t), intent(in), allocatable :: core_state
      real(default) :: alpha_qed
    end function prc_omega_get_alpha_qed
  end interface

contains

  subroutine omega_def_init (object, &
       model_name, prt_in, prt_out, &
       ovm, ufo, ufo_path, &
       restrictions, cms_scheme, &
       openmp_support, report_progress, write_phs_output, extra_options, diags, diags_color)
    class(omega_def_t), intent(out) :: object
    type(string_t), intent(in) :: model_name
    type(string_t), dimension(:), intent(in) :: prt_in
    type(string_t), dimension(:), intent(in) :: prt_out
    logical, intent(in) :: ovm
    logical, intent(in) :: ufo
    type(string_t), intent(in), optional :: ufo_path
    type(string_t), intent(in), optional :: restrictions
    logical, intent(in), optional :: cms_scheme
    logical, intent(in), optional :: openmp_support
    logical, intent(in), optional :: report_progress
    logical, intent(in), optional :: write_phs_output
    type(string_t), intent(in), optional :: extra_options
    logical, intent(in), optional :: diags, diags_color
    object%ufo = ufo
    object%ovm = ovm
    if (object%ufo) then
       if (object%ovm) then
          call msg_fatal ("Omega process: OVM method does not support UFO model")
       else
          allocate (omega_ufo_writer_t :: object%writer)
       end if
    else
       if (object%ovm) then
          allocate (omega_ovm_writer_t :: object%writer)
       else
          allocate (omega_omega_writer_t :: object%writer)
       end if
    end if
    select type (writer => object%writer)
    class is (omega_writer_t)
       call writer%init (model_name, prt_in, prt_out, &
            ufo_path, restrictions, cms_scheme, &
            openmp_support, report_progress, write_phs_output, extra_options, diags, diags_color)
    end select
  end subroutine omega_def_init

  subroutine omega_def_allocate_driver (object, driver, basename)
    class(omega_def_t), intent(in) :: object
    class(prc_core_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    allocate (omega_driver_t :: driver)
  end subroutine omega_def_allocate_driver

  subroutine omega_make_process_component (entry, component_index, &
         model_name, prt_in, prt_out, &
         ufo, ufo_path, restrictions, cms_scheme, &
         openmp_support, report_progress, write_omega_output, extra_options, diags, diags_color)
    class(process_def_entry_t), intent(inout) :: entry
    integer, intent(in) :: component_index
    type(string_t), intent(in) :: model_name
    type(string_t), dimension(:), intent(in) :: prt_in
    type(string_t), dimension(:), intent(in) :: prt_out
    logical, intent(in), optional :: ufo
    type(string_t), intent(in), optional :: ufo_path
    type(string_t), intent(in), optional :: restrictions
    logical, intent(in), optional :: cms_scheme
    logical, intent(in), optional :: openmp_support
    logical, intent(in), optional :: report_progress
    logical, intent(in), optional :: write_omega_output
    type(string_t), intent(in), optional :: extra_options
    logical, intent(in), optional :: diags, diags_color
    logical :: ufo_model
    class(prc_core_def_t), allocatable :: def
    ufo_model = .false.;  if (present (ufo))  ufo_model = ufo
    allocate (omega_def_t :: def)
    select type (def)
    class is (omega_def_t)
       call def%init (model_name, prt_in, prt_out, &
            .false., ufo_model, ufo_path, &
            restrictions, cms_scheme, &
            openmp_support, report_progress, write_omega_output, extra_options, diags, diags_color)
    end select
    call entry%import_component (component_index, &
         n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method = var_str ("omega"), &
         variant = def)
  end subroutine omega_make_process_component

  subroutine prc_omega_allocate_workspace (object, core_state)
    class(prc_omega_t), intent(in) :: object
    class(prc_core_state_t), intent(inout), allocatable :: core_state
    allocate (omega_state_t :: core_state)
  end subroutine prc_omega_allocate_workspace

  subroutine prc_omega_set_parameters (prc_omega, model, &
       helicity_selection, qcd, use_color_factors)
    class(prc_omega_t), intent(inout) :: prc_omega
    class(model_data_t), intent(in), target, optional :: model
    type(helicity_selection_t), intent(in), optional :: helicity_selection
    type(qcd_t), intent(in), optional :: qcd
    type(qed_t) :: qed
    logical, intent(in), optional :: use_color_factors
    if (present (model)) then
       if (.not. allocated (prc_omega%par)) &
            allocate (prc_omega%par (model%get_n_real ()))
       call model%real_parameters_to_array (prc_omega%par)
       prc_omega%scheme = model%get_scheme_num ()
       if (associated (model%get_par_data_ptr (var_str ('alpha_em_i')))) then
          allocate (alpha_qed_fixed_t :: qed%alpha)
          select type (alpha => qed%alpha)
          type is (alpha_qed_fixed_t)
             alpha%val = one / model%get_real (var_str ('alpha_em_i'))
          end select
       end if
       prc_omega%qed = qed
    end if
    if (present (helicity_selection)) then
       prc_omega%helicity_selection = helicity_selection
    end if
    if (present (qcd)) then
       prc_omega%qcd = qcd
    end if
    if (present (use_color_factors)) then
       prc_omega%use_color_factors = use_color_factors
    end if
  end subroutine prc_omega_set_parameters


end module prc_omega
