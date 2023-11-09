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
module prc_external

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds
  use iso_varying_string, string_t => varying_string
  use constants
  use os_interface
  use lorentz
  use interactions
  use sm_qcd
  use sm_qed
  use variables, only: var_list_t

  use model_data
  use prclib_interfaces
  use prc_core_def
  use prc_core

  use sf_base
  use sf_pdf_builtin, only: pdf_builtin_t
  use sf_lhapdf, only: lhapdf_t

  implicit none
  private

  public :: prc_external_state_t
  public :: prc_external_driver_t
  public :: prc_external_t
  public :: prc_external_def_t
  public :: prc_external_writer_t
  public :: prc_external_test_writer_t
  public :: prc_external_test_state_t
  public :: prc_external_test_driver_t
  public :: prc_external_test_def_t
  public :: prc_external_test_t

  type :: sf_handler_t
     integer :: initial_state_type = 0
     integer :: n_sf = -1
     real(default) :: val = one
  contains
    procedure :: init => sf_handler_init
    procedure :: init_dummy => sf_handler_init_dummy
    procedure :: apply_structure_functions => &
         sf_handler_apply_structure_functions
    procedure :: get_pdf => sf_handler_get_pdf
  end type sf_handler_t

  type, abstract, extends (prc_core_state_t) :: prc_external_state_t
    logical :: new_kinematics = .true.
    real(default) :: alpha_qcd = -1
    real(default) :: alpha_qed = -1
  contains
    procedure :: reset_new_kinematics => prc_external_state_reset_new_kinematics
  end type prc_external_state_t

  type, abstract, extends (prc_core_driver_t) :: prc_external_driver_t
     procedure(omega_update_alpha_s), nopass, pointer :: &
              update_alpha_s => null ()
     procedure(omega_is_allowed), nopass, pointer :: &
              is_allowed => null ()
  end type prc_external_driver_t

  type, abstract, extends (prc_core_t) :: prc_external_t
    type(qcd_t) :: qcd
    type(qed_t) :: qed
    integer :: n_flv = 1
    real(default), dimension(:), allocatable :: par
    integer :: scheme = 0
    type(sf_handler_t) :: sf_handler
    real(default) :: maximum_accuracy = 10000.0
  contains
    procedure, nopass :: needs_external_code => &
         prc_external_needs_external_code
    procedure :: get_n_flvs => prc_external_get_n_flvs
    procedure :: get_flv_state => prc_external_get_flv_state
    procedure :: compute_sqme => prc_external_compute_sqme
    procedure :: compute_sqme_virt => prc_external_compute_sqme_virt
    procedure :: compute_sqme_color_c => prc_external_compute_sqme_color_c
    procedure :: compute_alpha_s => prc_external_compute_alpha_s
    procedure :: get_alpha_s => prc_external_get_alpha_s
    procedure :: get_alpha_qed => prc_external_get_alpha_qed
    procedure :: is_allowed => prc_external_is_allowed
    procedure :: get_nflv => prc_external_get_nflv
    procedure :: compute_hard_kinematics => prc_external_compute_hard_kinematics
    procedure :: compute_eff_kinematics => prc_external_compute_eff_kinematics
    procedure :: set_parameters => prc_external_set_parameters
    procedure :: update_alpha_s => prc_external_update_alpha_s
    procedure :: init_sf_handler => prc_external_init_sf_handler
    procedure :: init_sf_handler_dummy => prc_external_init_sf_handler_dummy
    procedure :: apply_structure_functions => &
         prc_external_apply_structure_functions
    procedure :: get_sf_value => prc_external_get_sf_value
    procedure(prc_external_includes_polarization), deferred :: &
      includes_polarization
  end type prc_external_t

  type, abstract, extends (prc_core_def_t) :: prc_external_def_t
    type(string_t) :: basename
  contains
    procedure :: set_active_writer => prc_external_def_set_active_writer
    procedure, nopass :: get_features => prc_external_def_get_features
    procedure :: connect => prc_external_def_connect
    procedure :: omega_connect => prc_external_def_connect
    procedure, nopass :: needs_code => prc_external_def_needs_code
  end type prc_external_def_t

  type, abstract, extends (prc_writer_f_module_t) :: prc_external_writer_t
    type(string_t) :: model_name
    type(string_t) :: process_mode
    type(string_t) :: process_string
    type(string_t) :: restrictions
    integer :: n_in = 0
    integer :: n_out = 0
    logical :: active = .true.
    logical :: amp_triv = .true.
  contains
    procedure :: init => prc_external_writer_init
    procedure :: base_init => prc_external_writer_init
    procedure, nopass :: get_module_name => prc_external_writer_get_module_name
    procedure :: write_wrapper => prc_external_writer_write_wrapper
    procedure :: write_interface => prc_external_writer_write_interface
    procedure :: write_source_code => prc_external_writer_write_source_code
    procedure :: before_compile => prc_external_writer_before_compile
    procedure :: after_compile => prc_external_writer_after_compile
    procedure :: write_makefile_code => prc_external_writer_write_makefile_code
    procedure :: base_write_makefile_code => &
         prc_external_writer_write_makefile_code
    procedure, nopass:: get_procname => prc_external_writer_writer_get_procname
  end type prc_external_writer_t

  type, extends (prc_external_writer_t) :: prc_external_test_writer_t
  contains
    procedure, nopass :: type_name => prc_external_test_writer_type_name
  end type prc_external_test_writer_t

  type, extends (prc_external_state_t) :: prc_external_test_state_t
  contains
    procedure :: write => prc_external_test_state_write
  end type prc_external_test_state_t

  type, extends (prc_external_driver_t) :: prc_external_test_driver_t
  contains
    procedure, nopass :: type_name => prc_external_test_driver_type_name
  end type prc_external_test_driver_t

  type, extends (prc_external_def_t) :: prc_external_test_def_t
  contains
    procedure :: init => prc_external_test_def_init
    procedure, nopass :: type_string => prc_external_test_def_type_string
    procedure :: write => prc_external_test_def_write
    procedure :: read => prc_external_test_def_read
    procedure :: allocate_driver => prc_external_test_def_allocate_driver
  end type prc_external_test_def_t

  type, extends (prc_external_t) :: prc_external_test_t
  contains
    procedure :: write => prc_external_test_write
    procedure :: write_name => prc_external_test_write_name
    procedure :: compute_amplitude => prc_external_test_compute_amplitude
    procedure :: allocate_workspace => prc_external_test_allocate_workspace
    procedure :: includes_polarization => &
         prc_external_test_includes_polarization
    procedure :: prepare_external_code => &
         prc_external_test_prepare_external_code
  end type prc_external_test_t


  abstract interface
    function prc_external_includes_polarization (object) result (polarized)
      import
      logical :: polarized
      class(prc_external_t), intent(in) :: object
    end function prc_external_includes_polarization
  end interface

  abstract interface
     subroutine omega_update_alpha_s (alpha_s) bind(C)
       import
       real(c_default_float), intent(in) :: alpha_s
     end subroutine omega_update_alpha_s
  end interface

  abstract interface
     subroutine omega_is_allowed (flv, hel, col, flag) bind(C)
       import
       integer(c_int), intent(in) :: flv, hel, col
       logical(c_bool), intent(out) :: flag
     end subroutine omega_is_allowed
  end interface


  interface
    module subroutine sf_handler_init (sf_handler, sf_chain)
      class(sf_handler_t), intent(out) :: sf_handler
      type(sf_chain_instance_t), intent(in) :: sf_chain
    end subroutine sf_handler_init
    module subroutine sf_handler_init_dummy (sf_handler)
      class(sf_handler_t), intent(out) :: sf_handler
    end subroutine sf_handler_init_dummy
    module subroutine sf_handler_apply_structure_functions &
         (sf_handler, sf_chain, flavors)
       class(sf_handler_t), intent(inout) :: sf_handler
       type(sf_chain_instance_t), intent(in) :: sf_chain
       integer, intent(in), dimension(2) :: flavors
    end subroutine sf_handler_apply_structure_functions
    module function sf_handler_get_pdf &
         (sf_handler, sf_chain, i, flavor) result (f)
       real(default) :: f
       class(sf_handler_t), intent(in) :: sf_handler
       type(sf_chain_instance_t), intent(in) :: sf_chain
       integer, intent(in) :: i, flavor
    end function sf_handler_get_pdf
    module subroutine prc_external_state_reset_new_kinematics (object)
      class(prc_external_state_t), intent(inout) :: object
    end subroutine prc_external_state_reset_new_kinematics
    module function prc_external_needs_external_code () result (flag)
      logical :: flag
    end function prc_external_needs_external_code
    pure module function prc_external_get_n_flvs (object, i_flv) result (n)
      integer :: n
      class(prc_external_t), intent(in) :: object
      integer, intent(in) :: i_flv
    end function prc_external_get_n_flvs
    module function prc_external_get_flv_state (object, i_flv) result (flv)
      integer, dimension(:), allocatable :: flv
      class(prc_external_t), intent(in) :: object
      integer, intent(in) :: i_flv
    end function prc_external_get_flv_state
    module subroutine prc_external_compute_sqme (object, i_flv, i_hel, p, &
           ren_scale, sqme, bad_point)
       class(prc_external_t), intent(in) :: object
       integer, intent(in) :: i_flv, i_hel
       type(vector4_t), dimension(:), intent(in) :: p
       real(default), intent(in) :: ren_scale
       real(default), intent(out) :: sqme
       logical, intent(out) :: bad_point
    end subroutine prc_external_compute_sqme
    module subroutine prc_external_compute_sqme_virt (object, i_flv, i_hel, &
       p, ren_scale, es_scale, loop_method, sqme, bad_point)
      class(prc_external_t), intent(in) :: object
      integer, intent(in) :: i_flv, i_hel
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: ren_scale, es_scale
      integer, intent(in) :: loop_method
      logical, intent(out) :: bad_point
      real(default), dimension(4), intent(out) :: sqme
    end subroutine prc_external_compute_sqme_virt
    module subroutine prc_external_compute_sqme_color_c (object, i_flv, &
         i_hel, p, ren_scale, born_color_c, bad_point, born_out)
      class(prc_external_t), intent(inout) :: object
      integer, intent(in) :: i_flv, i_hel
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(in) :: ren_scale
      real(default), intent(inout), dimension(:,:) :: born_color_c
      logical, intent(out) :: bad_point
      real(default), intent(out), optional :: born_out
    end subroutine prc_external_compute_sqme_color_c
    module subroutine prc_external_compute_alpha_s &
         (object, core_state, ren_scale)
      class(prc_external_t), intent(in) :: object
      class(prc_external_state_t), intent(inout) :: core_state
      real(default), intent(in) :: ren_scale
    end subroutine prc_external_compute_alpha_s
    module function prc_external_get_alpha_s &
         (object, core_state) result (alpha_qcd)
      class(prc_external_t), intent(in) :: object
      class(prc_core_state_t), intent(in), allocatable :: core_state
      real(default) :: alpha_qcd
    end function prc_external_get_alpha_s
    module function prc_external_get_alpha_qed &
         (object, core_state) result (alpha_qed)
      class(prc_external_t), intent(in) :: object
      class(prc_core_state_t), intent(in), allocatable :: core_state
      real(default) :: alpha_qed
    end function prc_external_get_alpha_qed
    module function prc_external_is_allowed &
         (object, i_term, f, h, c) result (flag)
      class(prc_external_t), intent(in) :: object
      integer, intent(in) :: i_term, f, h, c
      logical :: flag
    end function prc_external_is_allowed
    module function prc_external_get_nflv (object) result (n_flv)
      class(prc_external_t), intent(in) :: object
      integer :: n_flv
    end function prc_external_get_nflv
    module subroutine prc_external_compute_hard_kinematics &
         (object, p_seed, i_term, int_hard, core_state)
      class(prc_external_t), intent(in) :: object
      type(vector4_t), dimension(:), intent(in) :: p_seed
      integer, intent(in) :: i_term
      type(interaction_t), intent(inout) :: int_hard
      class(prc_core_state_t), intent(inout), allocatable :: core_state
    end subroutine prc_external_compute_hard_kinematics
    module subroutine prc_external_compute_eff_kinematics &
         (object, i_term, int_hard, int_eff, core_state)
      class(prc_external_t), intent(in) :: object
      integer, intent(in) :: i_term
      type(interaction_t), intent(in) :: int_hard
      type(interaction_t), intent(inout) :: int_eff
      class(prc_core_state_t), intent(inout), allocatable :: core_state
    end subroutine prc_external_compute_eff_kinematics
    module subroutine prc_external_update_alpha_s (object, core_state, scale)
      class(prc_external_t), intent(in) :: object
      class(prc_core_state_t), intent(inout), allocatable :: core_state
      real(default), intent(in) :: scale
    end subroutine prc_external_update_alpha_s
    module subroutine prc_external_init_sf_handler (core, sf_chain)
       class(prc_external_t), intent(inout) :: core
       type(sf_chain_instance_t), intent(in) :: sf_chain
    end subroutine prc_external_init_sf_handler
    module subroutine prc_external_init_sf_handler_dummy (core)
       class(prc_external_t), intent(inout) :: core
    end subroutine prc_external_init_sf_handler_dummy
    module subroutine prc_external_apply_structure_functions &
         (core, sf_chain, flavors)
      class(prc_external_t), intent(inout) :: core
      type(sf_chain_instance_t), intent(in) :: sf_chain
      integer, dimension(2), intent(in) :: flavors
    end subroutine prc_external_apply_structure_functions
    module function prc_external_get_sf_value (core) result (val)
      real(default) :: val
      class(prc_external_t), intent(in) :: core
    end function prc_external_get_sf_value
    module subroutine prc_external_def_set_active_writer (def, active)
      class(prc_external_def_t), intent(inout) :: def
      logical, intent(in) :: active
    end subroutine prc_external_def_set_active_writer
    module subroutine prc_external_def_get_features (features)
      type(string_t), dimension(:), allocatable, intent(out) :: features
    end subroutine prc_external_def_get_features
    module subroutine prc_external_def_connect (def, lib_driver, i, proc_driver)
      class(prc_external_def_t), intent(in) :: def
      class(prclib_driver_t), intent(in) :: lib_driver
      integer, intent(in) :: i
      class(prc_core_driver_t), intent(inout) :: proc_driver
    end subroutine prc_external_def_connect
    module function prc_external_def_needs_code () result (flag)
      logical :: flag
    end function prc_external_def_needs_code
    pure module subroutine prc_external_writer_init &
         (writer, model_name, prt_in, prt_out, restrictions)
      class(prc_external_writer_t), intent(inout) :: writer
      type(string_t), intent(in) :: model_name
      type(string_t), dimension(:), intent(in) :: prt_in, prt_out
      type(string_t), intent(in), optional :: restrictions
    end subroutine prc_external_writer_init
    module function prc_external_writer_get_module_name (id) result (name)
      type(string_t) :: name
      type(string_t), intent(in) :: id
    end function prc_external_writer_get_module_name
    module subroutine prc_external_writer_write_wrapper &
         (writer, unit, id, feature)
      class(prc_external_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id, feature
    end subroutine prc_external_writer_write_wrapper
    module subroutine prc_external_writer_write_interface &
         (writer, unit, id, feature)
      class(prc_external_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
      type(string_t), intent(in) :: feature
    end subroutine prc_external_writer_write_interface
    module subroutine prc_external_writer_write_source_code (writer, id)
      class(prc_external_writer_t), intent(in) :: writer
      type(string_t), intent(in) :: id
    end subroutine prc_external_writer_write_source_code
    module subroutine prc_external_writer_before_compile (writer, id)
      class(prc_external_writer_t), intent(in) :: writer
      type(string_t), intent(in) :: id
    end subroutine prc_external_writer_before_compile
    module subroutine prc_external_writer_after_compile (writer, id)
      class(prc_external_writer_t), intent(in) :: writer
      type(string_t), intent(in) :: id
    end subroutine prc_external_writer_after_compile
    module subroutine prc_external_writer_write_makefile_code &
         (writer, unit, id, os_data, verbose, testflag)
      class(prc_external_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: verbose
      logical, intent(in), optional :: testflag
    end subroutine prc_external_writer_write_makefile_code
    module function prc_external_writer_writer_get_procname &
         (feature) result (name)
      type(string_t) :: name
      type(string_t), intent(in) :: feature
    end function prc_external_writer_writer_get_procname
    module function prc_external_test_writer_type_name () result (string)
      type(string_t) :: string
    end function prc_external_test_writer_type_name
    module subroutine prc_external_test_state_write (object, unit)
      class(prc_external_test_state_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_external_test_state_write
    module function prc_external_test_driver_type_name () result (type)
      type(string_t) :: type
    end function prc_external_test_driver_type_name
    module function prc_external_test_def_type_string () result (string)
      type(string_t) :: string
    end function prc_external_test_def_type_string
    module subroutine prc_external_test_def_write (object, unit)
      class(prc_external_test_def_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine prc_external_test_def_write
    module subroutine prc_external_test_def_read (object, unit)
      class(prc_external_test_def_t), intent(out) :: object
      integer, intent(in) :: unit
    end subroutine prc_external_test_def_read
    module subroutine prc_external_test_write (object, unit)
      class(prc_external_test_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_external_test_write
    module subroutine prc_external_test_write_name (object, unit)
      class(prc_external_test_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_external_test_write_name
    module function prc_external_test_compute_amplitude &
         (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
         core_state)  result (amp)
      class(prc_external_test_t), intent(in) :: object
      integer, intent(in) :: j
      type(vector4_t), dimension(:), intent(in) :: p
      integer, intent(in) :: f, h, c
      real(default), intent(in) :: fac_scale, ren_scale
      real(default), intent(in), allocatable :: alpha_qcd_forced
      class(prc_core_state_t), intent(inout), allocatable, optional :: &
           core_state
      complex(default) :: amp
    end function prc_external_test_compute_amplitude
    module function prc_external_test_includes_polarization &
         (object) result (polarized)
      logical :: polarized
      class(prc_external_test_t), intent(in) :: object
    end function prc_external_test_includes_polarization
    module subroutine prc_external_test_prepare_external_code &
         (core, flv_states, var_list, os_data, libname, model, i_core, is_nlo)
      class(prc_external_test_t), intent(inout) :: core
      integer, intent(in), dimension(:,:), allocatable :: flv_states
      type(var_list_t), intent(in) :: var_list
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in) :: libname
      type(model_data_t), intent(in), target :: model
      integer, intent(in) :: i_core
      logical, intent(in) :: is_nlo
    end subroutine prc_external_test_prepare_external_code
  end interface

contains

  subroutine prc_external_set_parameters (object, qcd, model)
    class(prc_external_t), intent(inout) :: object
    type(qcd_t), intent(in) :: qcd
    type(qed_t) :: qed
    class(model_data_t), intent(in), target, optional :: model
    object%qcd = qcd
    if (present (model)) then
       if (.not. allocated (object%par)) &
            allocate (object%par (model%get_n_real ()))
       call model%real_parameters_to_array (object%par)
       object%scheme = model%get_scheme_num ()
       if (associated (model%get_par_data_ptr (var_str ('alpha_em_i')))) then
          allocate (alpha_qed_fixed_t :: qed%alpha)
          select type (alpha => qed%alpha)
          type is (alpha_qed_fixed_t)
             alpha%val = one / model%get_real (var_str ('alpha_em_i'))
          end select
       end if
       object%qed = qed
    end if
  end subroutine prc_external_set_parameters

  subroutine prc_external_test_def_init (object, basename, model_name, &
       prt_in, prt_out)
    class(prc_external_test_def_t), intent(inout) :: object
    type(string_t), intent(in) :: basename, model_name
    type(string_t), dimension(:), intent(in) :: prt_in, prt_out
    object%basename = basename
    allocate (prc_external_test_writer_t :: object%writer)
    select type (writer => object%writer)
    type is (prc_external_test_writer_t)
       call writer%init (model_name, prt_in, prt_out)
    end select
  end subroutine prc_external_test_def_init

  subroutine prc_external_test_def_allocate_driver (object, driver, basename)
    class(prc_external_test_def_t), intent(in) :: object
    class(prc_core_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    if (.not. allocated (driver)) &
         allocate (prc_external_test_driver_t :: driver)
  end subroutine prc_external_test_def_allocate_driver

  subroutine prc_external_test_allocate_workspace (object, core_state)
    class(prc_external_test_t), intent(in) :: object
    class(prc_core_state_t), intent(inout), allocatable :: core_state
    allocate (prc_external_test_state_t :: core_state)
  end subroutine prc_external_test_allocate_workspace


end module prc_external
