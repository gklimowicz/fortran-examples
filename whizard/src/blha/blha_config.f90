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

module blha_config

  use kinds
  use iso_varying_string, string_t => varying_string
  use variables, only: var_list_t
  use model_data
  use beam_structures, only: beam_structure_t

  implicit none
  private

  public :: blha_configuration_t
  public :: blha_cfg_process_node_t
  public :: ew_scheme_string_to_int
  public :: correction_type_string_to_int
  public :: blha_flv_state_t
  public :: blha_master_t
  public :: blha_get_additional_suffix
  public :: blha_configuration_init
  public :: blha_configuration_append_processes
  public :: blha_configuration_set
  public :: blha_configuration_get_n_proc
  public :: blha_configuration_write

  integer, public, parameter :: &
       BLHA_CT_QCD = 1, BLHA_CT_EW = 2, BLHA_CT_OTHER = 3
  integer, public, parameter :: &
       BLHA_IRREG_CDR = 1, BLHA_IRREG_DRED = 2, BLHA_IRREG_THV = 3, &
       BLHA_IRREG_MREG = 4, BLHA_IRREG_OTHER = 5
  integer, public, parameter :: &
       BLHA_MPS_ONSHELL = 1, BLHA_MPS_OTHER = 2
  integer, public, parameter :: &
       BLHA_MODE_GOSAM = 1, BLHA_MODE_FEYNARTS = 2, BLHA_MODE_GENERIC = 3, &
       BLHA_MODE_OPENLOOPS = 4
  integer, public, parameter :: &
       BLHA_VERSION_1 = 1, BLHA_VERSION_2 = 2
  integer, public, parameter :: &
       BLHA_AMP_LOOP = 1, BLHA_AMP_COLOR_C = 2, BLHA_AMP_SPIN_C = 3, &
       BLHA_AMP_TREE = 4, BLHA_AMP_LOOPINDUCED = 5
  integer, public, parameter :: &
       BLHA_EW_INTERNAL = 0, &
       BLHA_EW_GF = 1, BLHA_EW_MZ = 2, BLHA_EW_MSBAR = 3, &
       BLHA_EW_0 = 4, BLHA_EW_RUN = 5
  integer, public, parameter :: &
       BLHA_WIDTH_COMPLEX = 1, BLHA_WIDTH_FIXED = 2, &
       BLHA_WIDTH_RUNNING = 3, BLHA_WIDTH_POLE = 4, &
       BLHA_WIDTH_DEFAULT = 5

  integer, parameter, public :: OLP_N_MASSIVE_PARTICLES = 12
  integer, dimension(OLP_N_MASSIVE_PARTICLES), public :: &
       OLP_MASSIVE_PARTICLES = [5, -5, 6, -6, 13, -13, 15, -15, 23, 24, -24, 25]
  integer, parameter :: OLP_HEL_UNPOLARIZED = 0

  integer, parameter :: N_KNOWN_SPECIAL_OL_METHODS = 3

  type :: blha_particle_string_element_t
     integer :: pdg = 0
     integer :: hel = OLP_HEL_UNPOLARIZED
     logical :: polarized = .false.
  contains
    generic :: init => init_default
    generic :: init => init_polarized
    procedure :: init_default => blha_particle_string_element_init_default
    procedure :: init_polarized => blha_particle_string_element_init_polarized
    generic :: write_pdg => write_pdg_unit
    generic :: write_pdg => write_pdg_character
    procedure :: write_pdg_unit => blha_particle_string_element_write_pdg_unit
    procedure :: write_pdg_character &
         => blha_particle_string_element_write_pdg_character
    generic :: write_helicity => write_helicity_unit
    generic :: write_helicity => write_helicity_character
    procedure :: write_helicity_unit &
         => blha_particle_string_element_write_helicity_unit
    procedure :: write_helicity_character &
         => blha_particle_string_element_write_helicity_character
  end type blha_particle_string_element_t

  type :: blha_cfg_process_node_t
     type(blha_particle_string_element_t), dimension(:), allocatable :: pdg_in, pdg_out
     integer, dimension(:), allocatable :: fingerprint
     integer :: nsub
     integer, dimension(:), allocatable :: ids
     integer :: amplitude_type
     type(blha_cfg_process_node_t), pointer :: next => null ()
  end type blha_cfg_process_node_t

  type :: blha_configuration_t
     type(string_t) :: name
     class(model_data_t), pointer :: model => null ()
     type(string_t) :: md5
     integer :: version = 2
     logical :: dirty = .false.
     integer :: n_proc = 0
     real(default) :: accuracy_target
     logical :: debug_unstable = .false.
     integer :: mode = BLHA_MODE_GENERIC
     logical :: polarized = .false.
     type(blha_cfg_process_node_t), pointer :: processes => null ()
     !integer, dimension(2) :: matrix_element_square_type = BLHA_MEST_SUM
     integer :: correction_type
     type(string_t) :: correction_type_other
     integer :: irreg = BLHA_IRREG_THV
     type(string_t) :: irreg_other
     integer :: massive_particle_scheme = BLHA_MPS_ONSHELL
     type(string_t) :: massive_particle_scheme_other
     type(string_t) :: model_file
     logical :: subdivide_subprocesses = .false.
     integer :: alphas_power = -1, alpha_power = -1
     integer :: ew_scheme = BLHA_EW_GF
     integer :: width_scheme = BLHA_WIDTH_DEFAULT
     logical :: openloops_use_cms = .false.
     integer :: openloops_phs_tolerance = 0
     type(string_t) :: openloops_extra_cmd
     type(string_t) :: openloops_allowed_libs
     integer :: openloops_stability_log = 0
     integer :: n_off_photons_is = 0
     integer :: n_off_photons_fs = 0
  end type blha_configuration_t

  type:: blha_flv_state_t
    integer, dimension(:), allocatable :: flavors
    integer :: flv_mult
    logical :: flv_real = .false.
  end type blha_flv_state_t

  type :: blha_master_t
    integer, dimension(5) :: blha_mode = BLHA_MODE_GENERIC
    logical :: compute_borns = .false.
    logical :: compute_real_trees = .false.
    logical :: compute_loops = .true.
    logical :: compute_correlations = .false.
    logical :: compute_dglap = .false.
    integer :: ew_scheme
    type(string_t), dimension(:), allocatable :: suffix
    type(blha_configuration_t), dimension(:), allocatable :: blha_cfg
    integer :: n_files = 0
    integer, dimension(:), allocatable :: i_file_to_nlo_index
  contains
    procedure :: set_methods => blha_master_set_methods
    procedure :: allocate_config_files => blha_master_allocate_config_files
    procedure :: set_ew_scheme => blha_master_set_ew_scheme
    procedure :: set_correction_type => blha_master_set_correction_type
    procedure :: set_photon_characteristics => blha_master_set_photon_characteristics
    procedure :: generate => blha_master_generate
    procedure :: generate_loop => blha_master_generate_loop
    procedure :: generate_correlation => blha_master_generate_correlation
    procedure :: generate_real_tree => blha_master_generate_real_tree
    procedure :: generate_born => blha_master_generate_born
    procedure :: generate_dglap => blha_master_generate_dglap
    procedure :: setup_additional_features => blha_master_setup_additional_features
    procedure :: set_gosam => blha_master_set_gosam
    procedure :: set_openloops => blha_master_set_openloops
    procedure :: set_polarization => blha_master_set_polarization
    procedure :: write_olp => blha_master_write_olp
    procedure :: final => blha_master_final
  end type blha_master_t






  interface
    module subroutine blha_particle_string_element_init_default (blha_p, id)
      class(blha_particle_string_element_t), intent(out) :: blha_p
      integer, intent(in) :: id
    end subroutine blha_particle_string_element_init_default
    module subroutine blha_particle_string_element_init_polarized (blha_p, id, hel)
      class(blha_particle_string_element_t), intent(out) :: blha_p
      integer, intent(in) :: id, hel
    end subroutine blha_particle_string_element_init_polarized
    module subroutine blha_particle_string_element_write_pdg_unit (blha_p, unit)
      class(blha_particle_string_element_t), intent(in) :: blha_p
      integer, intent(in), optional :: unit
    end subroutine blha_particle_string_element_write_pdg_unit
    module subroutine blha_particle_string_element_write_pdg_character (blha_p, c)
      class(blha_particle_string_element_t), intent(in) :: blha_p
      character(3), intent(inout) :: c
    end subroutine blha_particle_string_element_write_pdg_character
    module subroutine blha_particle_string_element_write_helicity_unit (blha_p, unit)
      class(blha_particle_string_element_t), intent(in) :: blha_p
      integer, intent(in), optional :: unit
    end subroutine blha_particle_string_element_write_helicity_unit
    module subroutine blha_particle_string_element_write_helicity_character (blha_p, c)
      class(blha_particle_string_element_t), intent(in) :: blha_p
      character(4), intent(inout) :: c
    end subroutine blha_particle_string_element_write_helicity_character
    module function ew_scheme_string_to_int (ew_scheme_str) result (ew_scheme_int)
      integer :: ew_scheme_int
      type(string_t), intent(in) :: ew_scheme_str
    end function ew_scheme_string_to_int
    module function correction_type_string_to_int &
         (correction_type_str) result (correction_type_int)
      integer :: correction_type_int
      type(string_t), intent(in) :: correction_type_str
    end function correction_type_string_to_int
    module subroutine blha_master_set_methods (master, is_nlo, var_list)
      class(blha_master_t), intent(inout) :: master
      logical, intent(in) :: is_nlo
      type(var_list_t), intent(in) :: var_list
    end subroutine blha_master_set_methods
    module subroutine blha_master_allocate_config_files (master)
      class(blha_master_t), intent(inout) :: master
    end subroutine blha_master_allocate_config_files
    module subroutine blha_master_set_ew_scheme (master, ew_scheme)
      class(blha_master_t), intent(inout) :: master
      type(string_t), intent(in) :: ew_scheme
    end subroutine blha_master_set_ew_scheme
    module subroutine blha_master_set_correction_type (master, correction_type_str)
      class(blha_master_t), intent(inout) :: master
      type(string_t), intent(in) :: correction_type_str
    end subroutine blha_master_set_correction_type
    module subroutine blha_master_set_photon_characteristics (master, flv_born, n_in)
      class(blha_master_t), intent(inout) :: master
      integer, dimension(:,:), intent(in) :: flv_born
      integer, intent(in) :: n_in
    end subroutine blha_master_set_photon_characteristics
    module subroutine blha_master_generate (master, basename, model, &
         n_in, alpha_power, alphas_power, flv_born, flv_real)
      class(blha_master_t), intent(inout) :: master
      type(string_t), intent(in) :: basename
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_in
      integer, intent(in) :: alpha_power, alphas_power
      integer, intent(in), dimension(:,:), allocatable :: flv_born, flv_real
    end subroutine blha_master_generate
    module subroutine blha_master_generate_loop (master, basename, model, n_in, &
         alpha_power, alphas_power, flv_born, i_file)
      class(blha_master_t), intent(inout) :: master
      type(string_t), intent(in) :: basename
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_in
      integer, intent(in) :: alpha_power, alphas_power
      integer, dimension(:,:), allocatable, intent(in) :: flv_born
      integer, intent(inout) :: i_file
    end subroutine blha_master_generate_loop
    module subroutine blha_master_generate_correlation (master, basename, model, n_in, &
         alpha_power, alphas_power, flv_born, i_file)
      class(blha_master_t), intent(inout) :: master
      type(string_t), intent(in) :: basename
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_in
      integer, intent(in) :: alpha_power, alphas_power
      integer, dimension(:,:), allocatable, intent(in) :: flv_born
      integer, intent(inout) :: i_file
    end subroutine blha_master_generate_correlation
    module subroutine blha_master_generate_real_tree (master, basename, model, n_in, &
         alpha_power, alphas_power, flv_real, i_file)
      class(blha_master_t), intent(inout) :: master
      type(string_t), intent(in) :: basename
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_in
      integer, intent(in) :: alpha_power, alphas_power
      integer, dimension(:,:), allocatable, intent(in) :: flv_real
      integer, intent(inout) :: i_file
    end subroutine blha_master_generate_real_tree
    module subroutine blha_master_generate_born (master, basename, model, n_in, &
         alpha_power, alphas_power, flv_born, i_file)
      class(blha_master_t), intent(inout) :: master
      type(string_t), intent(in) :: basename
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_in
      integer, intent(in) :: alpha_power, alphas_power
      integer, dimension(:,:), allocatable, intent(in) :: flv_born
      integer, intent(inout) :: i_file
    end subroutine blha_master_generate_born
    module subroutine blha_master_generate_dglap (master, basename, model, n_in, &
         alpha_power, alphas_power, flv_born, i_file)
      class(blha_master_t), intent(inout) :: master
      type(string_t), intent(in) :: basename
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_in
      integer, intent(in) :: alpha_power, alphas_power
      integer, dimension(:,:), allocatable, intent(in) :: flv_born
      integer, intent(inout) :: i_file
    end subroutine blha_master_generate_dglap
    module subroutine blha_master_setup_additional_features (master, &
           phs_tolerance, use_cms, stability_log, extra_cmd, &
           allowed_libs, beam_structure)
      class(blha_master_t), intent(inout) :: master
      integer, intent(in) :: phs_tolerance
      logical, intent(in) :: use_cms
      type(string_t), intent(in), optional :: extra_cmd, allowed_libs
      integer, intent(in) :: stability_log
      type(beam_structure_t), intent(in), optional :: beam_structure
    end subroutine blha_master_setup_additional_features
    module subroutine blha_master_set_gosam (master, i)
      class(blha_master_t), intent(inout) :: master
      integer, intent(in) :: i
    end subroutine blha_master_set_gosam
    module subroutine blha_master_set_openloops (master, i)
      class(blha_master_t), intent(inout) :: master
      integer, intent(in) :: i
    end subroutine blha_master_set_openloops
    module subroutine blha_master_set_polarization (master, i)
      class(blha_master_t), intent(inout) :: master
      integer, intent(in) :: i
    end subroutine blha_master_set_polarization
    module function blha_get_additional_suffix (base_suffix) result (suffix)
      type(string_t) :: suffix
      type(string_t), intent(in) :: base_suffix
    end function blha_get_additional_suffix
    module subroutine blha_master_write_olp (master, basename)
      class(blha_master_t), intent(in) :: master
      type(string_t), intent(in) :: basename
    end subroutine blha_master_write_olp
    module subroutine blha_master_final (master)
      class(blha_master_t), intent(inout) :: master
    end subroutine blha_master_final
    module subroutine blha_configuration_init (cfg, name, model, mode)
      type(blha_configuration_t), intent(inout) :: cfg
      type(string_t), intent(in) :: name
      class(model_data_t), target, intent(in) :: model
      integer, intent(in), optional :: mode
    end subroutine blha_configuration_init
    module subroutine blha_configuration_append_processes (cfg, n_in, flavor, amp_type)
      type(blha_configuration_t), intent(inout) :: cfg
      integer, intent(in) :: n_in
      type(blha_flv_state_t), dimension(:), intent(in) :: flavor
      integer, dimension(:), intent(in), optional :: amp_type
    end subroutine blha_configuration_append_processes
    module subroutine blha_configuration_set (cfg, &
         version, irreg, massive_particle_scheme, &
         model_file, alphas_power, alpha_power, ew_scheme, width_scheme, &
         accuracy, debug)
      type(blha_configuration_t), intent(inout) :: cfg
      integer, optional, intent(in) :: version
      integer, optional, intent(in) :: irreg
      integer, optional, intent(in) :: massive_particle_scheme
      type(string_t), optional, intent(in) :: model_file
      integer, optional, intent(in) :: alphas_power, alpha_power
      integer, optional, intent(in) :: ew_scheme
      integer, optional, intent(in) :: width_scheme
      real(default), optional, intent(in) :: accuracy
      logical, optional, intent(in) :: debug
    end subroutine blha_configuration_set
    module function blha_configuration_get_n_proc (cfg) result (n_proc)
      type(blha_configuration_t), intent(in) :: cfg
      integer :: n_proc
    end function blha_configuration_get_n_proc
    module subroutine blha_configuration_write (cfg, suffix, unit, internal, no_version)
      type(blha_configuration_t), intent(in) :: cfg
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: internal, no_version
      type(string_t), intent(in) :: suffix
    end subroutine blha_configuration_write
  end interface

end module blha_config

