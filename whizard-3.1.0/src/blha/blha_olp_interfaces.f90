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

module blha_olp_interfaces

  use, intrinsic :: iso_c_binding !NODEP!
  use, intrinsic :: iso_fortran_env

  use kinds
  use iso_varying_string, string_t => varying_string
  use os_interface
  use lorentz
  use interactions
  use model_data

  use prclib_interfaces
  use process_libraries
  use prc_core_def
  use prc_core
  use prc_external
  use blha_config

  implicit none
  private

  public :: blha_template_t
  public :: prc_blha_t
  public :: blha_driver_t
  public :: prc_blha_writer_t
  public :: blha_def_t
  public :: blha_state_t
  public :: olp_start
  public :: olp_eval
  public :: olp_info
  public :: olp_set_parameter
  public :: olp_eval2
  public :: olp_option
  public :: olp_polvec
  public :: olp_finalize
  public :: olp_print_parameter
  public :: blha_result_array_size
  public :: parameter_error_message
  public :: ew_parameter_error_message
  public :: blha_color_c_fill_diag
  public :: blha_color_c_fill_offdiag
  public :: blha_loop_positions

  integer, parameter, public :: OLP_PARAMETER_LIMIT = 10
  integer, parameter, public :: OLP_MOMENTUM_LIMIT = 50
  integer, parameter, public :: OLP_RESULTS_LIMIT = 60


  type :: blha_template_t
    integer :: I_BORN = 0
    integer :: I_REAL = 1
    integer :: I_LOOP = 2
    integer :: I_SUB = 3
    integer :: I_DGLAP = 4
    logical, dimension(0:4) :: compute_component
    logical :: include_polarizations = .false.
    logical :: switch_off_muon_yukawas = .false.
    logical :: use_internal_color_correlations = .true.
    real(default) :: external_top_yukawa = -1._default
    integer :: ew_scheme
    integer :: loop_method = BLHA_MODE_GENERIC
  contains
    procedure :: write => blha_template_write
    procedure :: get_n_hel => blha_template_get_n_hel
    procedure :: init => blha_template_init
    procedure :: set_born => blha_template_set_born
    procedure :: set_real_trees => blha_template_set_real_trees
    procedure :: set_loop => blha_template_set_loop
    procedure :: set_subtraction => blha_template_set_subtraction
    procedure :: set_dglap => blha_template_set_dglap
    procedure :: set_internal_color_correlations &
         => blha_template_set_internal_color_correlations
    procedure :: get_internal_color_correlations &
         => blha_template_get_internal_color_correlations
    procedure :: compute_born => blha_template_compute_born
    procedure :: compute_real_trees => blha_template_compute_real_trees
    procedure :: compute_loop => blha_template_compute_loop
    procedure :: compute_subtraction => blha_template_compute_subtraction
    procedure :: compute_dglap => blha_template_compute_dglap
    procedure :: set_loop_method => blha_template_set_loop_method
    procedure :: check => blha_template_check
    procedure :: reset => blha_template_reset
  end type blha_template_t

  type, abstract, extends (prc_external_t) :: prc_blha_t
    integer :: n_particles
    integer :: n_hel
    integer :: n_proc
    integer, dimension(:, :), allocatable :: i_tree, i_spin_c, i_color_c
    integer, dimension(:, :), allocatable :: i_virt
    integer, dimension(:, :), allocatable :: i_hel
    logical, dimension(5) :: ew_parameter_mask
    integer :: sqme_tree_pos
  contains
    procedure, nopass :: uses_blha => prc_blha_uses_blha
    procedure :: create_momentum_array => prc_blha_create_momentum_array
    procedure :: set_alpha_qed => prc_blha_set_alpha_qed
    procedure :: set_GF => prc_blha_set_GF
    procedure :: set_weinberg_angle => prc_blha_set_weinberg_angle
    procedure :: set_electroweak_parameters => &
         prc_blha_set_electroweak_parameters
    procedure :: read_contract_file => prc_blha_read_contract_file
    procedure :: print_parameter_file => prc_blha_print_parameter_file
    procedure :: compute_amplitude => prc_blha_compute_amplitude
    procedure :: init_blha => prc_blha_init_blha
    procedure :: set_mass_and_width => prc_blha_set_mass_and_width
    procedure :: set_particle_properties => prc_blha_set_particle_properties
    procedure :: init_ew_parameters => prc_blha_init_ew_parameters
    procedure :: compute_sqme_virt => prc_blha_compute_sqme_virt
    procedure :: compute_sqme => prc_blha_compute_sqme
    procedure :: compute_sqme_color_c_raw => prc_blha_compute_sqme_color_c_raw
    procedure :: compute_sqme_color_c => prc_blha_compute_sqme_color_c
    generic :: get_beam_helicities => get_beam_helicities_single
    generic :: get_beam_helicities => get_beam_helicities_array
    procedure :: get_beam_helicities_single => prc_blha_get_beam_helicities_single
    procedure :: get_beam_helicities_array => prc_blha_get_beam_helicities_array
    procedure :: includes_polarization => prc_blha_includes_polarization
    procedure :: set_equivalent_flv_hel_indices => prc_blha_set_equivalent_flv_hel_indices
    procedure(prc_blha_init_driver), deferred :: &
         init_driver
  end type prc_blha_t

  type, abstract, extends (prc_external_driver_t) :: blha_driver_t
    type(string_t) :: contract_file
    type(string_t) :: nlo_suffix
    logical :: include_polarizations = .false.
    logical :: switch_off_muon_yukawas = .false.
    real(default) :: external_top_yukawa = -1.0
    procedure(olp_start),nopass,  pointer :: &
         blha_olp_start => null ()
    procedure(olp_eval), nopass, pointer :: &
         blha_olp_eval => null()
    procedure(olp_info), nopass, pointer :: &
         blha_olp_info => null ()
    procedure(olp_set_parameter), nopass, pointer :: &
         blha_olp_set_parameter => null ()
    procedure(olp_eval2), nopass, pointer :: &
         blha_olp_eval2 => null ()
    procedure(olp_option), nopass, pointer :: &
         blha_olp_option => null ()
    procedure(olp_polvec), nopass, pointer :: &
         blha_olp_polvec => null ()
    procedure(olp_finalize), nopass, pointer :: &
         blha_olp_finalize => null ()
    procedure(olp_print_parameter), nopass, pointer :: &
         blha_olp_print_parameter => null ()
  contains
    procedure(blha_driver_set_GF), deferred :: &
         set_GF
    procedure(blha_driver_set_alpha_s), deferred :: &
         set_alpha_s
    procedure(blha_driver_set_weinberg_angle), deferred :: &
         set_weinberg_angle
    procedure(blha_driver_set_alpha_qed), deferred :: set_alpha_qed
    procedure(blha_driver_print_alpha_s), deferred :: &
         print_alpha_s
    procedure :: set_mass_and_width => blha_driver_set_mass_and_width
    procedure(blha_driver_init_dlaccess_to_library), deferred :: &
         init_dlaccess_to_library
    procedure :: load => blha_driver_load
    procedure :: read_contract_file => blha_driver_read_contract_file
  end type blha_driver_t

  type, abstract, extends (prc_external_writer_t) :: prc_blha_writer_t
    type(blha_configuration_t) :: blha_cfg
  contains
    procedure :: write => prc_blha_writer_write
    procedure :: get_process_string => prc_blha_writer_get_process_string
    procedure :: get_n_proc => prc_blha_writer_get_n_proc
  end type prc_blha_writer_t

  type, abstract, extends (prc_external_def_t) :: blha_def_t
    type(string_t) :: suffix
  contains
  
  end type blha_def_t

  type, abstract, extends (prc_external_state_t) :: blha_state_t
  contains
    procedure :: reset_new_kinematics => blha_state_reset_new_kinematics
  end type blha_state_t


  interface
    subroutine olp_start (contract_file_name, ierr) bind (C,name = "OLP_Start")
      import
      character(kind = c_char, len = 1), intent(in) :: contract_file_name
      integer(kind = c_int), intent(out) :: ierr
    end subroutine olp_start
  end interface

  interface
    subroutine olp_eval (label, momenta, mu, parameters, res) &
           bind (C, name = "OLP_EvalSubProcess")
      import
      integer(kind = c_int), value, intent(in) :: label
      real(kind = c_double), value, intent(in) :: mu
      real(kind = c_double), dimension(OLP_MOMENTUM_LIMIT), intent(in) :: &
           momenta
      real(kind = c_double), dimension(OLP_PARAMETER_LIMIT), intent(in) :: &
           parameters
      real(kind = c_double), dimension(OLP_RESULTS_LIMIT), intent(out) :: res
    end subroutine olp_eval
  end interface

  interface
    subroutine olp_info (olp_file, olp_version, message) bind(C)
      import
      character(kind = c_char), intent(inout), dimension(15) :: olp_file
      character(kind = c_char), intent(inout), dimension(15) :: olp_version
      character(kind = c_char), intent(inout), dimension(255) :: message
    end subroutine olp_info
  end interface

  interface
    subroutine olp_set_parameter &
           (variable_name, real_part, complex_part, success) bind(C)
      import
      character(kind = c_char,len = 1), intent(in) :: variable_name
      real(kind = c_double), intent(in) :: real_part, complex_part
      integer(kind = c_int), intent(out) :: success
    end subroutine olp_set_parameter
  end interface

  interface
    subroutine olp_eval2 (label, momenta, mu, res, acc) bind(C)
      import
      integer(kind = c_int), intent(in) :: label
      real(kind = c_double), intent(in) :: mu
      real(kind = c_double), dimension(OLP_MOMENTUM_LIMIT), intent(in) :: momenta
      real(kind = c_double), dimension(OLP_RESULTS_LIMIT), intent(out) :: res
      real(kind = c_double), intent(out) :: acc
    end subroutine olp_eval2
  end interface

  interface
    subroutine olp_option (line, stat) bind(C)
      import
      character(kind = c_char, len=1), intent(in) :: line
      integer(kind = c_int), intent(out) :: stat
    end subroutine
  end interface

  interface
    subroutine olp_polvec (p, q, eps) bind(C)
      import
      real(kind = c_double), dimension(0:3), intent(in) :: p, q
      real(kind = c_double), dimension(0:7), intent(out) :: eps
    end subroutine
  end interface

  interface
    subroutine olp_finalize () bind(C)
      import
    end subroutine olp_finalize
  end interface

  interface
    subroutine olp_print_parameter (filename) bind(C)
      import
      character(kind = c_char, len = 1), intent(in) :: filename
    end subroutine olp_print_parameter
  end interface

  abstract interface
    subroutine blha_driver_set_GF (driver, GF)
      import
      class(blha_driver_t), intent(inout) :: driver
      real(default), intent(in) :: GF
    end subroutine blha_driver_set_GF
  end interface

  abstract interface
    subroutine blha_driver_set_alpha_s (driver, alpha_s)
       import
       class(blha_driver_t), intent(in) :: driver
       real(default), intent(in) :: alpha_s
    end subroutine blha_driver_set_alpha_s
  end interface

  abstract interface
    subroutine blha_driver_set_weinberg_angle (driver, sw2)
      import
      class(blha_driver_t), intent(inout) :: driver
      real(default), intent(in) :: sw2
    end subroutine blha_driver_set_weinberg_angle
  end interface

  abstract interface
    subroutine blha_driver_set_alpha_qed (driver, alpha)
      import
      class(blha_driver_t), intent(inout) :: driver
      real(default), intent(in) :: alpha
    end subroutine blha_driver_set_alpha_qed
  end interface

  abstract interface
    subroutine blha_driver_print_alpha_s (object)
      import
      class(blha_driver_t), intent(in) :: object
    end subroutine blha_driver_print_alpha_s
  end interface

  abstract interface
    subroutine blha_driver_init_dlaccess_to_library &
           (object, os_data, dlaccess, success)
      import
      class(blha_driver_t), intent(in) :: object
      type(os_data_t), intent(in) :: os_data
      type(dlaccess_t), intent(out) :: dlaccess
      logical, intent(out) :: success
    end subroutine blha_driver_init_dlaccess_to_library
  end interface

  abstract interface
    subroutine prc_blha_init_driver (object, os_data)
      import
      class(prc_blha_t), intent(inout) :: object
      type(os_data_t), intent(in) :: os_data
    end subroutine prc_blha_init_driver
  end interface


  interface
    module subroutine blha_template_write (blha_template, unit)
      class(blha_template_t), intent(in) :: blha_template
      integer, intent(in), optional :: unit
    end subroutine blha_template_write
    module function blha_template_get_n_hel (blha_template, pdg, model) result (n_hel)
      class(blha_template_t), intent(in) :: blha_template
      integer, dimension(:), intent(in) :: pdg
      class(model_data_t), intent(in), target :: model
      integer :: n_hel
    end function blha_template_get_n_hel
    module function prc_blha_uses_blha () result (flag)
      logical :: flag
    end function prc_blha_uses_blha
    module subroutine blha_state_reset_new_kinematics (object)
      class(blha_state_t), intent(inout) :: object
    end subroutine blha_state_reset_new_kinematics
    pure module function blha_result_array_size &
         (n_part, amp_type) result (rsize)
      integer, intent(in) :: n_part, amp_type
      integer :: rsize
    end function blha_result_array_size
    module function prc_blha_create_momentum_array (object, p) result (mom)
      class(prc_blha_t), intent(in) :: object
      type(vector4_t), intent(in), dimension(:) :: p
      real(double), dimension(5*object%n_particles) :: mom
    end function prc_blha_create_momentum_array
    module subroutine blha_template_init (template, requires_polarizations, &
           switch_off_muon_yukawas, external_top_yukawa, ew_scheme)
      class(blha_template_t), intent(inout) :: template
      logical, intent(in) :: requires_polarizations, switch_off_muon_yukawas
      real(default), intent(in) :: external_top_yukawa
      type(string_t), intent(in) :: ew_scheme
    end subroutine blha_template_init
    module subroutine blha_template_set_born (template)
      class(blha_template_t), intent(inout) :: template
    end subroutine blha_template_set_born
    module subroutine blha_template_set_real_trees (template)
      class(blha_template_t), intent(inout) :: template
    end subroutine blha_template_set_real_trees
    module subroutine blha_template_set_loop (template)
      class(blha_template_t), intent(inout) :: template
    end subroutine blha_template_set_loop
    module subroutine blha_template_set_subtraction (template)
      class(blha_template_t), intent(inout) :: template
    end subroutine blha_template_set_subtraction
    module subroutine blha_template_set_dglap (template)
      class(blha_template_t), intent(inout) :: template
    end subroutine blha_template_set_dglap
    module subroutine blha_template_set_internal_color_correlations (template)
      class(blha_template_t), intent(inout) :: template
    end subroutine blha_template_set_internal_color_correlations
    pure module function blha_template_get_internal_color_correlations &
         (template) result (val)
      logical :: val
      class(blha_template_t), intent(in) :: template
    end function blha_template_get_internal_color_correlations
    pure module function blha_template_compute_born (template) result (val)
      class(blha_template_t), intent(in) :: template
      logical :: val
    end function blha_template_compute_born
    pure module function blha_template_compute_real_trees (template) result (val)
      class(blha_template_t), intent(in) :: template
      logical :: val
    end function blha_template_compute_real_trees
    pure module function blha_template_compute_loop (template) result (val)
      class(blha_template_t), intent(in) :: template
      logical :: val
    end function blha_template_compute_loop
    pure module function blha_template_compute_subtraction (template) result (val)
      class(blha_template_t), intent(in) :: template
      logical :: val
    end function blha_template_compute_subtraction
    pure module function blha_template_compute_dglap (template) result (val)
      class(blha_template_t), intent(in) :: template
      logical :: val
    end function blha_template_compute_dglap
    module subroutine blha_template_set_loop_method (template, master)
      class(blha_template_t), intent(inout) :: template
      class(blha_master_t), intent(in) :: master
    end subroutine blha_template_set_loop_method
    module function blha_template_check (template) result (val)
      class(blha_template_t), intent(in) :: template
      logical :: val
    end function blha_template_check
    module subroutine blha_template_reset (template)
      class(blha_template_t), intent(inout) :: template
    end subroutine blha_template_reset
    module subroutine prc_blha_writer_write (writer, unit)
      class(prc_blha_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
    end subroutine prc_blha_writer_write
    module function prc_blha_writer_get_process_string (writer) result (s_proc)
      class(prc_blha_writer_t), intent(in) :: writer
      type(string_t) :: s_proc
    end function prc_blha_writer_get_process_string
    module function prc_blha_writer_get_n_proc (writer) result (n_proc)
      class(prc_blha_writer_t), intent(in) :: writer
      integer :: n_proc
    end function prc_blha_writer_get_n_proc
    module subroutine parameter_error_message (par, subr)
      type(string_t), intent(in) :: par, subr
    end subroutine parameter_error_message
    module subroutine ew_parameter_error_message (par)
      type(string_t), intent(in) :: par
    end subroutine ew_parameter_error_message
    module subroutine blha_driver_set_mass_and_width &
         (driver, i_pdg, mass, width)
      class(blha_driver_t), intent(inout) :: driver
      integer, intent(in) :: i_pdg
      real(default), intent(in), optional :: mass
      real(default), intent(in), optional :: width
    end subroutine blha_driver_set_mass_and_width
    module subroutine blha_driver_load (object, os_data, success)
      class(blha_driver_t), intent(inout) :: object
      type(os_data_t), intent(in) :: os_data
      logical, intent(out) :: success
    end subroutine blha_driver_load
    module subroutine blha_driver_read_contract_file (driver, flavors, &
         amp_type, flv_index, hel_index, label, helicities)
      class(blha_driver_t), intent(inout) :: driver
      integer, intent(in), dimension(:,:) :: flavors
      integer, intent(out), dimension(:), allocatable :: amp_type, &
           flv_index, hel_index, label
      integer, intent(out), dimension(:,:) :: helicities
    end subroutine blha_driver_read_contract_file
    module subroutine prc_blha_set_alpha_qed (object, model)
      class(prc_blha_t), intent(inout) :: object
      type(model_data_t), intent(in), target :: model
    end subroutine prc_blha_set_alpha_qed
    module subroutine prc_blha_set_GF (object, model)
      class(prc_blha_t), intent(inout) :: object
      type(model_data_t), intent(in), target :: model
    end subroutine prc_blha_set_GF
    module subroutine prc_blha_set_weinberg_angle (object, model)
      class(prc_blha_t), intent(inout) :: object
      type(model_data_t), intent(in), target :: model
    end subroutine prc_blha_set_weinberg_angle
    module subroutine prc_blha_set_electroweak_parameters (object, model)
      class(prc_blha_t), intent(inout) :: object
      type(model_data_t), intent(in), target :: model
    end subroutine prc_blha_set_electroweak_parameters
    module subroutine prc_blha_read_contract_file (object, flavors)
      class(prc_blha_t), intent(inout) :: object
      integer, intent(in), dimension(:,:) :: flavors
    end subroutine prc_blha_read_contract_file
    module subroutine prc_blha_print_parameter_file (object, i_component)
      class(prc_blha_t), intent(in) :: object
      integer, intent(in) :: i_component
    end subroutine prc_blha_print_parameter_file
    module function prc_blha_compute_amplitude &
         (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
         core_state)  result (amp)
      class(prc_blha_t), intent(in) :: object
      integer, intent(in) :: j
      type(vector4_t), dimension(:), intent(in) :: p
      integer, intent(in) :: f, h, c
      real(default), intent(in) :: fac_scale, ren_scale
      real(default), intent(in), allocatable :: alpha_qcd_forced
      class(prc_core_state_t), intent(inout), allocatable, optional :: core_state
      complex(default) :: amp
    end function prc_blha_compute_amplitude
    module subroutine prc_blha_init_blha (object, blha_template, n_in, &
         n_particles, n_flv, n_hel)
      class(prc_blha_t), intent(inout) :: object
      type(blha_template_t), intent(in) :: blha_template
      integer, intent(in) :: n_in, n_particles, n_flv, n_hel
    end subroutine prc_blha_init_blha
    module subroutine prc_blha_set_mass_and_width (object, i_pdg, mass, width)
      class(prc_blha_t), intent(inout) :: object
      integer, intent(in) :: i_pdg
      real(default), intent(in) :: mass, width
    end subroutine prc_blha_set_mass_and_width
    module subroutine prc_blha_set_particle_properties (object, model)
      class(prc_blha_t), intent(inout) :: object
      class(model_data_t), intent(in), target :: model
    end subroutine prc_blha_set_particle_properties
    module subroutine prc_blha_init_ew_parameters (object, ew_scheme)
      class(prc_blha_t), intent(inout) :: object
      integer, intent(in) :: ew_scheme
    end subroutine prc_blha_init_ew_parameters
    module subroutine prc_blha_compute_sqme_virt (object, &
         i_flv, i_hel, p, ren_scale, es_scale, loop_method, sqme, bad_point)
      class(prc_blha_t), intent(in) :: object
      integer, intent(in) :: i_flv, i_hel
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: ren_scale, es_scale
      integer, intent(in) :: loop_method
      real(default), dimension(4), intent(out) :: sqme
      logical, intent(out) :: bad_point
    end subroutine prc_blha_compute_sqme_virt
    module subroutine prc_blha_compute_sqme (object, i_flv, i_hel, p, &
         ren_scale, sqme, bad_point)
      class(prc_blha_t), intent(in) :: object
      integer, intent(in) :: i_flv, i_hel
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(in) :: ren_scale
      real(default), intent(out) :: sqme
      logical, intent(out) :: bad_point
    end subroutine prc_blha_compute_sqme
    module subroutine blha_color_c_fill_diag &
         (sqme_born, flavors, sqme_color_c, special_case)
      real(default), intent(in) :: sqme_born
      integer, intent(in), dimension(:) :: flavors
      logical, intent(in), optional :: special_case
      real(default), intent(inout), dimension(:,:) :: sqme_color_c
    end subroutine blha_color_c_fill_diag
    module subroutine blha_color_c_fill_offdiag &
         (n, r, sqme_color_c, offset, n_flv)
      integer, intent(in) :: n
      real(default), intent(in), dimension(:) :: r
      real(default), intent(inout), dimension(:,:) :: sqme_color_c
      integer, intent(in), optional :: offset, n_flv
    end subroutine blha_color_c_fill_offdiag
    module subroutine prc_blha_compute_sqme_color_c_raw &
         (object, i_flv, i_hel, p, ren_scale, rr, bad_point)
      class(prc_blha_t), intent(in) :: object
      integer, intent(in) :: i_flv, i_hel
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(in) :: ren_scale
      real(default), intent(out), dimension(:) :: rr
      logical, intent(out) :: bad_point
    end subroutine prc_blha_compute_sqme_color_c_raw
    module subroutine prc_blha_compute_sqme_color_c &
         (object, i_flv, i_hel, p, ren_scale, born_color_c, bad_point, born_out)
      class(prc_blha_t), intent(inout) :: object
      integer, intent(in) :: i_flv, i_hel
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(in) :: ren_scale
      real(default), intent(inout), dimension(:,:) :: born_color_c
      real(default), intent(out), optional :: born_out
      logical, intent(out) :: bad_point
    end subroutine prc_blha_compute_sqme_color_c
    module function prc_blha_get_beam_helicities_single &
         (object, i, invert_second) result (hel)
      integer, dimension(:), allocatable :: hel
      class(prc_blha_t), intent(in) :: object
      logical, intent(in), optional :: invert_second
      integer, intent(in) :: i
    end function prc_blha_get_beam_helicities_single
    module function prc_blha_get_beam_helicities_array &
         (object, invert_second) result (hel)
      integer, dimension(:,:), allocatable :: hel
      class(prc_blha_t), intent(in) :: object
      logical, intent(in), optional :: invert_second
    end function prc_blha_get_beam_helicities_array
    module function prc_blha_includes_polarization (object) result (polarized)
      logical :: polarized
      class(prc_blha_t), intent(in) :: object
    end function prc_blha_includes_polarization
    module subroutine prc_blha_set_equivalent_flv_hel_indices (object)
      class(prc_blha_t), intent(inout) :: object
    end subroutine prc_blha_set_equivalent_flv_hel_indices
    recursive module function blha_loop_positions (i_flv, n_sub) result (index)
      integer :: index
      integer, intent(in) :: i_flv, n_sub
    end function blha_loop_positions
  end interface

end module blha_olp_interfaces

