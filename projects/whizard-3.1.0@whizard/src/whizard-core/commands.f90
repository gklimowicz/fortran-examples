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

module commands

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use diagnostics
  use lexers
  use syntax_rules
  use parser
  use variables, only: var_list_t, V_NONE, V_LOG, V_INT, V_REAL, V_CMPLX, V_STR, V_PDG
  use eval_trees
  use polarizations
  use rt_data

  implicit none
  private

  public :: get_prclib_static
  public :: command_list_t
  public :: syntax_cmd_list
  public :: syntax_cmd_list_init
  public :: syntax_cmd_list_final
  public :: syntax_cmd_list_write
  public :: lexer_init_cmd_list

  type, abstract :: command_t
     type(parse_node_t), pointer :: pn => null ()
     class(command_t), pointer :: next => null ()
     type(parse_node_t), pointer :: pn_opt => null ()
     type(command_list_t), pointer :: options => null ()
     type(rt_data_t), pointer :: local => null ()
   contains
     procedure :: final => command_final
     procedure (command_write), deferred :: write
     procedure (command_compile), deferred :: compile
     procedure (command_execute), deferred :: execute
     procedure :: write_options => command_write_options
     procedure :: compile_options => command_compile_options
     procedure :: execute_options => cmd_execute_options
     procedure :: reset_options => cmd_reset_options
  end type command_t

  type, extends (command_t) :: cmd_model_t
     private
     type(string_t) :: name
     type(string_t) :: scheme
     logical :: ufo_model = .false.
     logical :: ufo_path_set = .false.
     type(string_t) :: ufo_path
   contains
     procedure :: write => cmd_model_write
     procedure :: compile => cmd_model_compile
     procedure :: execute => cmd_model_execute
  end type cmd_model_t

  type, extends (command_t) :: cmd_library_t
     private
     type(string_t) :: name
   contains
     procedure :: write => cmd_library_write
     procedure :: compile => cmd_library_compile
     procedure :: execute => cmd_library_execute
  end type cmd_library_t

  type, extends (command_t) :: cmd_process_t
     private
     type(string_t) :: id
     integer :: n_in  = 0
     type(parse_node_p), dimension(:), allocatable :: pn_pdg_in
     type(parse_node_t), pointer :: pn_out => null ()
   contains
     procedure :: write => cmd_process_write
     procedure :: compile => cmd_process_compile
     procedure :: execute => cmd_process_execute
  end type cmd_process_t

  type, extends (command_t) :: cmd_nlo_t
    private
    integer, dimension(:), allocatable :: nlo_component
  contains
      procedure :: write => cmd_nlo_write
      procedure :: compile => cmd_nlo_compile
      procedure :: execute => cmd_nlo_execute
  end type cmd_nlo_t

  type, extends (command_t) :: cmd_compile_t
     private
     type(string_t), dimension(:), allocatable :: libname
     logical :: make_executable = .false.
     type(string_t) :: exec_name
   contains
     procedure :: write => cmd_compile_write
     procedure :: compile => cmd_compile_compile
     procedure :: execute => cmd_compile_execute
  end type cmd_compile_t

  type, extends (command_t) :: cmd_exec_t
     private
     type(parse_node_t), pointer :: pn_command => null ()
   contains
     procedure :: write => cmd_exec_write
     procedure :: compile => cmd_exec_compile
     procedure :: execute => cmd_exec_execute
  end type cmd_exec_t

  type, extends (command_t) :: cmd_var_t
     private
     type(string_t) :: name
     integer :: type = V_NONE
     type(parse_node_t), pointer :: pn_value => null ()
     logical :: is_intrinsic = .false.
     logical :: is_model_var = .false.
   contains
     procedure :: write => cmd_var_write
     procedure :: compile => cmd_var_compile
     procedure :: execute => cmd_var_execute
     procedure :: set_value => cmd_var_set_value
  end type cmd_var_t

  type, extends (command_t) :: cmd_slha_t
     private
     type(string_t) :: file
     logical :: write_mode = .false.
   contains
     procedure :: write => cmd_slha_write
     procedure :: compile => cmd_slha_compile
     procedure :: execute => cmd_slha_execute
  end type cmd_slha_t

  type, extends (command_t) :: cmd_show_t
     private
     type(string_t), dimension(:), allocatable :: name
   contains
     procedure :: write => cmd_show_write
     procedure :: compile => cmd_show_compile
     procedure :: execute => cmd_show_execute
  end type cmd_show_t

  type, extends (command_t) :: cmd_clear_t
     private
     type(string_t), dimension(:), allocatable :: name
   contains
     procedure :: write => cmd_clear_write
     procedure :: compile => cmd_clear_compile
     procedure :: execute => cmd_clear_execute
  end type cmd_clear_t

  type, extends (command_t) :: cmd_expect_t
     private
     type(parse_node_t), pointer :: pn_lexpr => null ()
   contains
     procedure :: write => cmd_expect_write
     procedure :: compile => cmd_expect_compile
     procedure :: execute => cmd_expect_execute
  end type cmd_expect_t

  type, extends (command_t) :: cmd_beams_t
     private
     integer :: n_in = 0
     type(parse_node_p), dimension(:), allocatable :: pn_pdg
     integer :: n_sf_record = 0
     integer, dimension(:), allocatable :: n_entry
     type(parse_node_p), dimension(:,:), allocatable :: pn_sf_entry
   contains
     procedure :: write => cmd_beams_write
     procedure :: compile => cmd_beams_compile
     procedure :: execute => cmd_beams_execute
  end type cmd_beams_t

  type :: sentry_expr_t
     type(parse_node_p), dimension(:), allocatable :: expr
   contains
     procedure :: compile => sentry_expr_compile
     procedure :: evaluate => sentry_expr_evaluate
  end type sentry_expr_t

  type :: smatrix_expr_t
     type(sentry_expr_t), dimension(:), allocatable :: entry
   contains
     procedure :: compile => smatrix_expr_compile
     procedure :: evaluate => smatrix_expr_evaluate
  end type smatrix_expr_t

  type, extends (command_t) :: cmd_beams_pol_density_t
     private
     integer :: n_in = 0
     type(smatrix_expr_t), dimension(:), allocatable :: smatrix
   contains
     procedure :: write => cmd_beams_pol_density_write
     procedure :: compile => cmd_beams_pol_density_compile
     procedure :: execute => cmd_beams_pol_density_execute
  end type cmd_beams_pol_density_t

  type, extends (command_t) :: cmd_beams_pol_fraction_t
     private
     integer :: n_in = 0
     type(parse_node_p), dimension(:), allocatable :: expr
   contains
     procedure :: write => cmd_beams_pol_fraction_write
     procedure :: compile => cmd_beams_pol_fraction_compile
     procedure :: execute => cmd_beams_pol_fraction_execute
  end type cmd_beams_pol_fraction_t

  type, extends (cmd_beams_pol_fraction_t) :: cmd_beams_momentum_t
   contains
     procedure :: write => cmd_beams_momentum_write
     procedure :: execute => cmd_beams_momentum_execute
  end type cmd_beams_momentum_t

  type, extends (cmd_beams_pol_fraction_t) :: cmd_beams_theta_t
   contains
     procedure :: write => cmd_beams_theta_write
     procedure :: execute => cmd_beams_theta_execute
  end type cmd_beams_theta_t

  type, extends (cmd_beams_pol_fraction_t) :: cmd_beams_phi_t
   contains
     procedure :: write => cmd_beams_phi_write
     procedure :: execute => cmd_beams_phi_execute
  end type cmd_beams_phi_t

  type, extends (command_t) :: cmd_cuts_t
     private
     type(parse_node_t), pointer :: pn_lexpr => null ()
   contains
     procedure :: write => cmd_cuts_write
     procedure :: compile => cmd_cuts_compile
     procedure :: execute => cmd_cuts_execute
  end type cmd_cuts_t

  type, extends (command_t) :: cmd_scale_t
     private
     type(parse_node_t), pointer :: pn_expr => null ()
   contains
     procedure :: write => cmd_scale_write
     procedure :: compile => cmd_scale_compile
     procedure :: execute => cmd_scale_execute
  end type cmd_scale_t

  type, extends (command_t) :: cmd_fac_scale_t
     private
     type(parse_node_t), pointer :: pn_expr => null ()
   contains
     procedure :: write => cmd_fac_scale_write
     procedure :: compile => cmd_fac_scale_compile
     procedure :: execute => cmd_fac_scale_execute
  end type cmd_fac_scale_t

  type, extends (command_t) :: cmd_ren_scale_t
     private
     type(parse_node_t), pointer :: pn_expr => null ()
   contains
     procedure :: write => cmd_ren_scale_write
     procedure :: compile => cmd_ren_scale_compile
     procedure :: execute => cmd_ren_scale_execute
  end type cmd_ren_scale_t

  type, extends (command_t) :: cmd_weight_t
     private
     type(parse_node_t), pointer :: pn_expr => null ()
   contains
     procedure :: write => cmd_weight_write
     procedure :: compile => cmd_weight_compile
     procedure :: execute => cmd_weight_execute
  end type cmd_weight_t

  type, extends (command_t) :: cmd_selection_t
     private
     type(parse_node_t), pointer :: pn_expr => null ()
   contains
     procedure :: write => cmd_selection_write
     procedure :: compile => cmd_selection_compile
     procedure :: execute => cmd_selection_execute
  end type cmd_selection_t

  type, extends (command_t) :: cmd_reweight_t
     private
     type(parse_node_t), pointer :: pn_expr => null ()
   contains
     procedure :: write => cmd_reweight_write
     procedure :: compile => cmd_reweight_compile
     procedure :: execute => cmd_reweight_execute
  end type cmd_reweight_t

  type, extends (command_t) :: cmd_alt_setup_t
     private
     type(parse_node_p), dimension(:), allocatable :: setup
   contains
     procedure :: write => cmd_alt_setup_write
     procedure :: compile => cmd_alt_setup_compile
     procedure :: execute => cmd_alt_setup_execute
  end type cmd_alt_setup_t

  type, extends (command_t) :: cmd_integrate_t
     private
     integer :: n_proc = 0
     type(string_t), dimension(:), allocatable :: process_id
   contains
     procedure :: write => cmd_integrate_write
     procedure :: compile => cmd_integrate_compile
     procedure :: execute => cmd_integrate_execute
  end type cmd_integrate_t

  type, extends (command_t) :: cmd_observable_t
     private
     type(string_t) :: id
   contains
     procedure :: write => cmd_observable_write
     procedure :: compile => cmd_observable_compile
     procedure :: execute => cmd_observable_execute
  end type cmd_observable_t

  type, extends (command_t) :: cmd_histogram_t
     private
     type(string_t) :: id
     type(parse_node_t), pointer :: pn_lower_bound => null ()
     type(parse_node_t), pointer :: pn_upper_bound => null ()
     type(parse_node_t), pointer :: pn_bin_width => null ()
   contains
     procedure :: write => cmd_histogram_write
     procedure :: compile => cmd_histogram_compile
     procedure :: execute => cmd_histogram_execute
  end type cmd_histogram_t

  type, extends (command_t) :: cmd_plot_t
     private
     type(string_t) :: id
   contains
     procedure :: write => cmd_plot_write
     procedure :: compile => cmd_plot_compile
     procedure :: init => cmd_plot_init
     procedure :: execute => cmd_plot_execute
  end type cmd_plot_t

  type, extends (command_t) :: cmd_graph_t
     private
     type(string_t) :: id
     integer :: n_elements = 0
     type(cmd_plot_t), dimension(:), allocatable :: el
     type(string_t), dimension(:), allocatable :: element_id
   contains
     procedure :: write => cmd_graph_write
     procedure :: compile => cmd_graph_compile
     procedure :: execute => cmd_graph_execute
  end type cmd_graph_t

  type :: analysis_id_t
    type(string_t) :: tag
    type(parse_node_t), pointer :: pn_sexpr => null ()
  end type analysis_id_t

  type, extends (command_t) :: cmd_analysis_t
     private
     type(parse_node_t), pointer :: pn_lexpr => null ()
   contains
     procedure :: write => cmd_analysis_write
     procedure :: compile => cmd_analysis_compile
     procedure :: execute => cmd_analysis_execute
  end type cmd_analysis_t

  type, extends (command_t) :: cmd_write_analysis_t
     private
     type(analysis_id_t), dimension(:), allocatable :: id
     type(string_t), dimension(:), allocatable :: tag
   contains
     procedure :: write => cmd_write_analysis_write
     procedure :: compile => cmd_write_analysis_compile
     procedure :: execute => cmd_write_analysis_execute
  end type cmd_write_analysis_t

  type, extends (command_t) :: cmd_compile_analysis_t
     private
     type(analysis_id_t), dimension(:), allocatable :: id
     type(string_t), dimension(:), allocatable :: tag
   contains
     procedure :: write => cmd_compile_analysis_write
     procedure :: compile => cmd_compile_analysis_compile
     procedure :: execute => cmd_compile_analysis_execute
  end type cmd_compile_analysis_t

  type, extends (command_t) :: cmd_open_out_t
     private
     type(parse_node_t), pointer :: file_expr => null ()
   contains
     procedure :: write => cmd_open_out_write
     procedure :: compile => cmd_open_out_compile
     procedure :: execute => cmd_open_out_execute
  end type cmd_open_out_t

  type, extends (cmd_open_out_t) :: cmd_close_out_t
     private
   contains
     procedure :: execute => cmd_close_out_execute
  end type cmd_close_out_t

  type, extends (command_t) :: cmd_printf_t
     private
     type(parse_node_t), pointer :: sexpr => null ()
     type(parse_node_t), pointer :: sprintf_fun => null ()
     type(parse_node_t), pointer :: sprintf_clause => null ()
     type(parse_node_t), pointer :: sprintf => null ()
   contains
     procedure :: final => cmd_printf_final
     procedure :: write => cmd_printf_write
     procedure :: compile => cmd_printf_compile
     procedure :: execute => cmd_printf_execute
  end type cmd_printf_t

  type, extends (command_t) :: cmd_record_t
     private
     type(parse_node_t), pointer :: pn_lexpr => null ()
   contains
     procedure :: write => cmd_record_write
     procedure :: compile => cmd_record_compile
     procedure :: execute => cmd_record_execute
  end type cmd_record_t

  type, extends (command_t) :: cmd_unstable_t
     private
     integer :: n_proc = 0
     type(string_t), dimension(:), allocatable :: process_id
     type(parse_node_t), pointer :: pn_prt_in => null ()
   contains
     procedure :: write => cmd_unstable_write
     procedure :: compile => cmd_unstable_compile
     procedure :: execute => cmd_unstable_execute
  end type cmd_unstable_t

  type, extends (command_t) :: cmd_stable_t
     private
     type(parse_node_p), dimension(:), allocatable :: pn_pdg
   contains
     procedure :: write => cmd_stable_write
     procedure :: compile => cmd_stable_compile
     procedure :: execute => cmd_stable_execute
  end type cmd_stable_t

  type, extends (cmd_stable_t) :: cmd_polarized_t
   contains
     procedure :: write => cmd_polarized_write
     procedure :: execute => cmd_polarized_execute
  end type cmd_polarized_t

  type, extends (cmd_stable_t) :: cmd_unpolarized_t
   contains
     procedure :: write => cmd_unpolarized_write
     procedure :: execute => cmd_unpolarized_execute
  end type cmd_unpolarized_t

  type, extends (command_t) :: cmd_sample_format_t
     private
     type(string_t), dimension(:), allocatable :: format
   contains
     procedure :: write => cmd_sample_format_write
     procedure :: compile => cmd_sample_format_compile
     procedure :: execute => cmd_sample_format_execute
  end type cmd_sample_format_t

  type, extends (command_t) :: cmd_simulate_t
     ! not private anymore as required by the whizard-c-interface
     integer :: n_proc = 0
     type(string_t), dimension(:), allocatable :: process_id
   contains
     procedure :: write => cmd_simulate_write
     procedure :: compile => cmd_simulate_compile
     procedure :: execute => cmd_simulate_execute
  end type cmd_simulate_t

  type, extends (command_t) :: cmd_rescan_t
     ! private
     type(parse_node_t), pointer :: pn_filename => null ()
     integer :: n_proc = 0
     type(string_t), dimension(:), allocatable :: process_id
   contains
     procedure :: write => cmd_rescan_write
     procedure :: compile => cmd_rescan_compile
     procedure :: execute => cmd_rescan_execute
  end type cmd_rescan_t

  type, extends (command_t) :: cmd_iterations_t
     private
     integer :: n_pass = 0
     type(parse_node_p), dimension(:), allocatable :: pn_expr_n_it
     type(parse_node_p), dimension(:), allocatable :: pn_expr_n_calls
     type(parse_node_p), dimension(:), allocatable :: pn_sexpr_adapt
   contains
     procedure :: write => cmd_iterations_write
     procedure :: compile => cmd_iterations_compile
     procedure :: execute => cmd_iterations_execute
  end type cmd_iterations_t

  type, abstract :: range_t
     type(parse_node_t), pointer :: pn_expr => null ()
     type(parse_node_t), pointer :: pn_term => null ()
     type(parse_node_t), pointer :: pn_factor => null ()
     type(parse_node_t), pointer :: pn_value => null ()
     type(parse_node_t), pointer :: pn_literal => null ()
     type(parse_node_t), pointer :: pn_beg => null ()
     type(parse_node_t), pointer :: pn_end => null ()
     type(parse_node_t), pointer :: pn_step => null ()
     type(eval_tree_t) :: expr_beg
     type(eval_tree_t) :: expr_end
     type(eval_tree_t) :: expr_step
     integer :: step_mode = 0
     integer :: n_step = 0
   contains
     procedure :: final => range_final
     procedure (range_write), deferred :: write
     procedure :: base_write => range_write
     procedure :: init => range_init
     procedure :: create_value_node => range_create_value_node
     procedure :: compile => range_compile
     procedure (range_evaluate), deferred :: evaluate
     procedure :: get_n_iterations => range_get_n_iterations
     procedure (range_set_value), deferred :: set_value
  end type range_t

  type, extends (range_t) :: range_int_t
     integer :: i_beg = 0
     integer :: i_end = 0
     integer :: i_step = 0
   contains
     procedure :: write => range_int_write
     procedure :: evaluate => range_int_evaluate
     procedure :: set_value => range_int_set_value
  end type range_int_t

  type, extends (range_t) :: range_real_t
     real(default) :: r_beg = 0
     real(default) :: r_end = 0
     real(default) :: r_step = 0
     real(default) :: lr_beg  = 0
     real(default) :: lr_end  = 0
     real(default) :: lr_step = 0
   contains
     procedure :: write => range_real_write
     procedure :: evaluate => range_real_evaluate
     procedure :: set_value => range_real_set_value
  end type range_real_t

  type, extends (command_t) :: cmd_scan_t
     private
     type(string_t) :: name
     integer :: n_values = 0
     type(parse_node_p), dimension(:), allocatable :: scan_cmd
     class(range_t), dimension(:), allocatable :: range
   contains
     procedure :: final => cmd_scan_final
     procedure :: write => cmd_scan_write
     procedure :: compile => cmd_scan_compile
     procedure :: execute => cmd_scan_execute
  end type cmd_scan_t

  type, extends (command_t) :: cmd_if_t
     private
     type(parse_node_t), pointer :: pn_if_lexpr => null ()
     type(command_list_t), pointer :: if_body => null ()
     type(cmd_if_t), dimension(:), pointer :: elsif_cmd => null ()
     type(command_list_t), pointer :: else_body => null ()
   contains
     procedure :: final => cmd_if_final
     procedure :: write => cmd_if_write
     procedure :: compile => cmd_if_compile
     procedure :: execute => cmd_if_execute
  end type cmd_if_t

  type, extends (command_t) :: cmd_include_t
     private
     type(string_t) :: file
     type(command_list_t), pointer :: command_list => null ()
     type(parse_tree_t) :: parse_tree
   contains
     procedure :: final => cmd_include_final
     procedure :: write => cmd_include_write
     procedure :: compile => cmd_include_compile
     procedure :: execute => cmd_include_execute
  end type cmd_include_t

  type, extends (command_t) :: cmd_export_t
     private
     type(string_t), dimension(:), allocatable :: name
   contains
     procedure :: write => cmd_export_write
     procedure :: compile => cmd_export_compile
     procedure :: execute => cmd_export_execute
  end type cmd_export_t

  type, extends (command_t) :: cmd_quit_t
     private
     logical :: has_code = .false.
     type(parse_node_t), pointer :: pn_code_expr => null ()
   contains
     procedure :: write => cmd_quit_write
     procedure :: compile => cmd_quit_compile
     procedure :: execute => cmd_quit_execute
  end type cmd_quit_t

  type :: command_list_t
     ! not private anymore as required by the whizard-c-interface
     class(command_t), pointer :: first => null ()
     class(command_t), pointer :: last => null ()
   contains
     procedure :: write => command_list_write
     procedure :: append => command_list_append
     procedure :: final => command_list_final
     procedure :: compile => command_list_compile
     procedure :: execute => command_list_execute
  end type command_list_t


  type(syntax_t), target, save :: syntax_cmd_list


  integer, parameter, public :: SHOW_BUFFER_SIZE = 4096
  character(*), parameter, public :: &
       DEFAULT_ANALYSIS_FILENAME = "whizard_analysis.dat"
  character(len=1), dimension(2), parameter, public :: &
       FORBIDDEN_ENDINGS1 = [ "o", "a" ]
  character(len=2), dimension(6), parameter, public :: &
       FORBIDDEN_ENDINGS2 = [ "mp", "ps", "vg", "pg", "lo", "la" ]
  character(len=3), dimension(20), parameter, public :: &
       FORBIDDEN_ENDINGS3 = [ "aux", "dvi", "evt", "evx", "f03", "f90", &
          "f95", "log", "ltp", "mod", "mpx", "olc", "olp", "pdf", "phs", &
          "sin", "sub", "tex", "vg2", "vgx" ]

  integer, parameter :: STEP_NONE = 0
  integer, parameter :: STEP_ADD = 1
  integer, parameter :: STEP_SUB = 2
  integer, parameter :: STEP_MUL = 3
  integer, parameter :: STEP_DIV = 4
  integer, parameter :: STEP_COMP_ADD = 11
  integer, parameter :: STEP_COMP_MUL = 13

  abstract interface
     subroutine command_write (cmd, unit, indent)
       import
       class(command_t), intent(in) :: cmd
       integer, intent(in), optional :: unit, indent
     end subroutine command_write
  end interface

  abstract interface
     subroutine command_compile (cmd, global)
       import
       class(command_t), intent(inout) :: cmd
       type(rt_data_t), intent(inout), target :: global
     end subroutine command_compile
  end interface

  abstract interface
     subroutine command_execute (cmd, global)
       import
       class(command_t), intent(inout) :: cmd
       type(rt_data_t), intent(inout), target :: global
     end subroutine command_execute
  end interface

  interface
     subroutine get_prclib_static (libname)
       import
       type(string_t), dimension(:), intent(inout), allocatable :: libname
     end subroutine get_prclib_static
  end interface

  abstract interface
     subroutine range_evaluate (range)
       import
       class(range_t), intent(inout) :: range
     end subroutine range_evaluate
  end interface

  abstract interface
     subroutine range_set_value (range, i)
       import
       class(range_t), intent(inout) :: range
       integer, intent(in) :: i
     end subroutine range_set_value
  end interface


  interface
    recursive module subroutine command_final (cmd)
      class(command_t), intent(inout) :: cmd
    end subroutine command_final
    recursive module subroutine command_write_options (cmd, unit, indent)
      class(command_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine command_write_options
    recursive module subroutine command_compile_options (cmd, global)
      class(command_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine command_compile_options
    recursive module subroutine cmd_execute_options (cmd, global)
      class(command_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_execute_options
    module subroutine cmd_reset_options (cmd, global)
      class(command_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_reset_options
    module subroutine cmd_model_write (cmd, unit, indent)
      class(cmd_model_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_model_write
    module subroutine cmd_model_compile (cmd, global)
      class(cmd_model_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_model_compile
    module subroutine cmd_model_execute (cmd, global)
      class(cmd_model_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_model_execute
    module subroutine cmd_library_write (cmd, unit, indent)
      class(cmd_library_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_library_write
    module subroutine cmd_library_compile (cmd, global)
      class(cmd_library_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_library_compile
    module subroutine cmd_library_execute (cmd, global)
      class(cmd_library_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_library_execute
    module subroutine cmd_process_write (cmd, unit, indent)
      class(cmd_process_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_process_write
    module subroutine cmd_process_compile (cmd, global)
      class(cmd_process_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_process_compile
    module subroutine cmd_process_execute (cmd, global)
      class(cmd_process_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_process_execute
    module subroutine cmd_nlo_write (cmd, unit, indent)
      class(cmd_nlo_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_nlo_write
    module subroutine cmd_nlo_compile (cmd, global)
      class(cmd_nlo_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_nlo_compile
    module subroutine cmd_nlo_execute (cmd, global)
      class(cmd_nlo_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_nlo_execute
    module subroutine cmd_compile_write (cmd, unit, indent)
      class(cmd_compile_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_compile_write
    module subroutine cmd_compile_compile (cmd, global)
      class(cmd_compile_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_compile_compile
    module subroutine cmd_compile_execute (cmd, global)
      class(cmd_compile_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_compile_execute
    module subroutine cmd_exec_write (cmd, unit, indent)
      class(cmd_exec_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_exec_write
    module subroutine cmd_exec_compile (cmd, global)
      class(cmd_exec_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_exec_compile
    module subroutine cmd_exec_execute (cmd, global)
      class(cmd_exec_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_exec_execute
    module subroutine cmd_var_write (cmd, unit, indent)
      class(cmd_var_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_var_write
    module subroutine cmd_var_compile (cmd, global)
      class(cmd_var_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_var_compile
    module subroutine cmd_var_execute (cmd, global)
      class(cmd_var_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_var_execute
    module subroutine cmd_var_set_value (var, var_list, verbose, model_name)
      class(cmd_var_t), intent(inout) :: var
      type(var_list_t), intent(inout), target :: var_list
      logical, intent(in), optional :: verbose
      type(string_t), intent(in), optional :: model_name
    end subroutine cmd_var_set_value
    module subroutine cmd_slha_write (cmd, unit, indent)
      class(cmd_slha_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_slha_write
    module subroutine cmd_slha_compile (cmd, global)
      class(cmd_slha_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_slha_compile
    module subroutine cmd_slha_execute (cmd, global)
      class(cmd_slha_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_slha_execute
    module subroutine cmd_show_write (cmd, unit, indent)
      class(cmd_show_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_show_write
    module subroutine cmd_show_compile (cmd, global)
      class(cmd_show_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_show_compile
    module subroutine cmd_show_execute (cmd, global)
      class(cmd_show_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_show_execute
    module subroutine cmd_clear_write (cmd, unit, indent)
      class(cmd_clear_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_clear_write
    module subroutine cmd_clear_compile (cmd, global)
      class(cmd_clear_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_clear_compile
    module subroutine cmd_clear_execute (cmd, global)
      class(cmd_clear_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_clear_execute
    module subroutine cmd_expect_write (cmd, unit, indent)
      class(cmd_expect_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_expect_write
    module subroutine cmd_expect_compile (cmd, global)
      class(cmd_expect_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_expect_compile
    module subroutine cmd_expect_execute (cmd, global)
      class(cmd_expect_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_expect_execute
    module subroutine cmd_beams_write (cmd, unit, indent)
      class(cmd_beams_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_beams_write
    module subroutine cmd_beams_compile (cmd, global)
      class(cmd_beams_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_beams_compile
    module subroutine cmd_beams_execute (cmd, global)
      class(cmd_beams_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_beams_execute
    module subroutine sentry_expr_compile (sentry, pn)
      class(sentry_expr_t), intent(out) :: sentry
      type(parse_node_t), intent(in), target :: pn
    end subroutine sentry_expr_compile
    module subroutine sentry_expr_evaluate (sentry, index, value, global)
      class(sentry_expr_t), intent(inout) :: sentry
      integer, dimension(:), intent(out) :: index
      complex(default), intent(out) :: value
      type(rt_data_t), intent(in), target :: global
    end subroutine sentry_expr_evaluate
    module subroutine smatrix_expr_compile (smatrix_expr, pn)
      class(smatrix_expr_t), intent(out) :: smatrix_expr
      type(parse_node_t), intent(in), target :: pn
    end subroutine smatrix_expr_compile
    module subroutine smatrix_expr_evaluate (smatrix_expr, smatrix, global)
      class(smatrix_expr_t), intent(inout) :: smatrix_expr
      type(smatrix_t), intent(out) :: smatrix
      type(rt_data_t), intent(in), target :: global
    end subroutine smatrix_expr_evaluate
    module subroutine cmd_beams_pol_density_write (cmd, unit, indent)
      class(cmd_beams_pol_density_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_beams_pol_density_write
    module subroutine cmd_beams_pol_density_compile (cmd, global)
      class(cmd_beams_pol_density_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_beams_pol_density_compile
    module subroutine cmd_beams_pol_density_execute (cmd, global)
      class(cmd_beams_pol_density_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_beams_pol_density_execute
    module subroutine cmd_beams_pol_fraction_write (cmd, unit, indent)
      class(cmd_beams_pol_fraction_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_beams_pol_fraction_write
    module subroutine cmd_beams_pol_fraction_compile (cmd, global)
      class(cmd_beams_pol_fraction_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_beams_pol_fraction_compile
    module subroutine cmd_beams_pol_fraction_execute (cmd, global)
      class(cmd_beams_pol_fraction_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_beams_pol_fraction_execute
    module subroutine cmd_beams_momentum_write (cmd, unit, indent)
      class(cmd_beams_momentum_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_beams_momentum_write
    module subroutine cmd_beams_momentum_execute (cmd, global)
      class(cmd_beams_momentum_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_beams_momentum_execute
    module subroutine cmd_beams_theta_write (cmd, unit, indent)
      class(cmd_beams_theta_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_beams_theta_write
    module subroutine cmd_beams_phi_write (cmd, unit, indent)
      class(cmd_beams_phi_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_beams_phi_write
    module subroutine cmd_beams_theta_execute (cmd, global)
      class(cmd_beams_theta_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_beams_theta_execute
    module subroutine cmd_beams_phi_execute (cmd, global)
      class(cmd_beams_phi_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_beams_phi_execute
    module subroutine cmd_cuts_write (cmd, unit, indent)
      class(cmd_cuts_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_cuts_write
    module subroutine cmd_cuts_compile (cmd, global)
      class(cmd_cuts_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_cuts_compile
    module subroutine cmd_cuts_execute (cmd, global)
      class(cmd_cuts_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_cuts_execute
    module subroutine cmd_scale_write (cmd, unit, indent)
      class(cmd_scale_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_scale_write
    module subroutine cmd_fac_scale_write (cmd, unit, indent)
      class(cmd_fac_scale_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_fac_scale_write
    module subroutine cmd_ren_scale_write (cmd, unit, indent)
      class(cmd_ren_scale_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_ren_scale_write
    module subroutine cmd_scale_compile (cmd, global)
      class(cmd_scale_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_scale_compile
    module subroutine cmd_fac_scale_compile (cmd, global)
      class(cmd_fac_scale_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_fac_scale_compile
    module subroutine cmd_ren_scale_compile (cmd, global)
      class(cmd_ren_scale_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_ren_scale_compile
    module subroutine cmd_scale_execute (cmd, global)
      class(cmd_scale_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_scale_execute
    module subroutine cmd_fac_scale_execute (cmd, global)
      class(cmd_fac_scale_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_fac_scale_execute
    module subroutine cmd_ren_scale_execute (cmd, global)
      class(cmd_ren_scale_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_ren_scale_execute
    module subroutine cmd_weight_write (cmd, unit, indent)
      class(cmd_weight_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_weight_write
    module subroutine cmd_weight_compile (cmd, global)
      class(cmd_weight_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_weight_compile
    module subroutine cmd_weight_execute (cmd, global)
      class(cmd_weight_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_weight_execute
    module subroutine cmd_selection_write (cmd, unit, indent)
      class(cmd_selection_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_selection_write
    module subroutine cmd_selection_compile (cmd, global)
      class(cmd_selection_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_selection_compile
    module subroutine cmd_selection_execute (cmd, global)
      class(cmd_selection_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_selection_execute
    module subroutine cmd_reweight_write (cmd, unit, indent)
      class(cmd_reweight_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_reweight_write
    module subroutine cmd_reweight_compile (cmd, global)
      class(cmd_reweight_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_reweight_compile
    module subroutine cmd_reweight_execute (cmd, global)
      class(cmd_reweight_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_reweight_execute
    module subroutine cmd_alt_setup_write (cmd, unit, indent)
      class(cmd_alt_setup_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_alt_setup_write
    module subroutine cmd_alt_setup_compile (cmd, global)
      class(cmd_alt_setup_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_alt_setup_compile
    module subroutine cmd_alt_setup_execute (cmd, global)
      class(cmd_alt_setup_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_alt_setup_execute
    module subroutine cmd_integrate_write (cmd, unit, indent)
      class(cmd_integrate_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_integrate_write
    module subroutine cmd_integrate_compile (cmd, global)
      class(cmd_integrate_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_integrate_compile
    module subroutine cmd_integrate_execute (cmd, global)
      class(cmd_integrate_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_integrate_execute
    module subroutine cmd_observable_write (cmd, unit, indent)
      class(cmd_observable_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_observable_write
    module subroutine cmd_observable_compile (cmd, global)
      class(cmd_observable_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_observable_compile
    module subroutine cmd_observable_execute (cmd, global)
      class(cmd_observable_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_observable_execute
    module subroutine cmd_histogram_write (cmd, unit, indent)
      class(cmd_histogram_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_histogram_write
    module subroutine cmd_histogram_compile (cmd, global)
      class(cmd_histogram_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_histogram_compile
    module subroutine cmd_histogram_execute (cmd, global)
      class(cmd_histogram_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_histogram_execute
    module subroutine cmd_plot_write (cmd, unit, indent)
      class(cmd_plot_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_plot_write
    module subroutine cmd_plot_compile (cmd, global)
      class(cmd_plot_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_plot_compile
    module subroutine cmd_plot_init (plot, pn_tag, global)
      class(cmd_plot_t), intent(inout) :: plot
      type(parse_node_t), intent(in), pointer :: pn_tag
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_plot_init
    module subroutine cmd_plot_execute (cmd, global)
      class(cmd_plot_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_plot_execute
    module subroutine cmd_graph_write (cmd, unit, indent)
      class(cmd_graph_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_graph_write
    module subroutine cmd_graph_compile (cmd, global)
      class(cmd_graph_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_graph_compile
    module subroutine cmd_graph_execute (cmd, global)
      class(cmd_graph_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_graph_execute
    module subroutine cmd_analysis_write (cmd, unit, indent)
      class(cmd_analysis_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_analysis_write
    module subroutine cmd_analysis_compile (cmd, global)
      class(cmd_analysis_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_analysis_compile
    module subroutine cmd_analysis_execute (cmd, global)
      class(cmd_analysis_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_analysis_execute
    module subroutine cmd_write_analysis_write (cmd, unit, indent)
      class(cmd_write_analysis_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_write_analysis_write
    module subroutine cmd_write_analysis_compile (cmd, global)
      class(cmd_write_analysis_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_write_analysis_compile
    module subroutine cmd_write_analysis_execute (cmd, global)
      class(cmd_write_analysis_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_write_analysis_execute
    module subroutine cmd_compile_analysis_write (cmd, unit, indent)
      class(cmd_compile_analysis_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_compile_analysis_write
    module subroutine cmd_compile_analysis_compile (cmd, global)
      class(cmd_compile_analysis_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_compile_analysis_compile
    module subroutine cmd_compile_analysis_execute (cmd, global)
      class(cmd_compile_analysis_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_compile_analysis_execute
    module subroutine cmd_open_out_write (cmd, unit, indent)
      class(cmd_open_out_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_open_out_write
    module subroutine cmd_open_out_compile (cmd, global)
      class(cmd_open_out_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_open_out_compile
    module subroutine cmd_open_out_execute (cmd, global)
      class(cmd_open_out_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_open_out_execute
    module subroutine cmd_close_out_execute (cmd, global)
      class(cmd_close_out_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_close_out_execute
    module subroutine cmd_printf_final (cmd)
      class(cmd_printf_t), intent(inout) :: cmd
    end subroutine cmd_printf_final
    module subroutine cmd_printf_write (cmd, unit, indent)
      class(cmd_printf_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_printf_write
    module subroutine cmd_printf_compile (cmd, global)
      class(cmd_printf_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_printf_compile
    module subroutine cmd_printf_execute (cmd, global)
      class(cmd_printf_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_printf_execute
    module subroutine cmd_record_write (cmd, unit, indent)
      class(cmd_record_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_record_write
    module subroutine cmd_record_compile (cmd, global)
      class(cmd_record_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_record_compile
    module subroutine cmd_record_execute (cmd, global)
      class(cmd_record_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_record_execute
    module subroutine cmd_unstable_write (cmd, unit, indent)
      class(cmd_unstable_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_unstable_write
    module subroutine cmd_unstable_compile (cmd, global)
      class(cmd_unstable_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_unstable_compile
    module subroutine cmd_unstable_execute (cmd, global)
      class(cmd_unstable_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_unstable_execute
    module subroutine cmd_stable_write (cmd, unit, indent)
      class(cmd_stable_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_stable_write
    module subroutine cmd_stable_compile (cmd, global)
      class(cmd_stable_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_stable_compile
    module subroutine cmd_stable_execute (cmd, global)
      class(cmd_stable_t), intent(inout) :: cmd
      type(rt_data_t), target, intent(inout) :: global
    end subroutine cmd_stable_execute
    module subroutine cmd_polarized_write (cmd, unit, indent)
      class(cmd_polarized_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_polarized_write
    module subroutine cmd_unpolarized_write (cmd, unit, indent)
      class(cmd_unpolarized_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_unpolarized_write
    module subroutine cmd_polarized_execute (cmd, global)
      class(cmd_polarized_t), intent(inout) :: cmd
      type(rt_data_t), target, intent(inout) :: global
    end subroutine cmd_polarized_execute
    module subroutine cmd_unpolarized_execute (cmd, global)
      class(cmd_unpolarized_t), intent(inout) :: cmd
      type(rt_data_t), target, intent(inout) :: global
    end subroutine cmd_unpolarized_execute
    module subroutine cmd_sample_format_write (cmd, unit, indent)
      class(cmd_sample_format_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_sample_format_write
    module subroutine cmd_sample_format_compile (cmd, global)
      class(cmd_sample_format_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_sample_format_compile
    module subroutine cmd_sample_format_execute (cmd, global)
      class(cmd_sample_format_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_sample_format_execute
    module subroutine cmd_simulate_write (cmd, unit, indent)
      class(cmd_simulate_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_simulate_write
    module subroutine cmd_simulate_compile (cmd, global)
      class(cmd_simulate_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_simulate_compile
    module subroutine cmd_simulate_execute (cmd, global)
      class(cmd_simulate_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_simulate_execute
    module subroutine cmd_rescan_write (cmd, unit, indent)
      class(cmd_rescan_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_rescan_write
    module subroutine cmd_rescan_compile (cmd, global)
      class(cmd_rescan_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_rescan_compile
    module subroutine cmd_rescan_execute (cmd, global)
      class(cmd_rescan_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_rescan_execute
    module subroutine cmd_iterations_write (cmd, unit, indent)
      class(cmd_iterations_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_iterations_write
    module subroutine cmd_iterations_compile (cmd, global)
      class(cmd_iterations_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_iterations_compile
  module subroutine cmd_iterations_execute (cmd, global)
    class(cmd_iterations_t), intent(inout) :: cmd
    type(rt_data_t), intent(inout), target :: global
  end subroutine cmd_iterations_execute
    module subroutine range_final (object)
      class(range_t), intent(inout) :: object
    end subroutine range_final
    module subroutine range_write (object, unit)
      class(range_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine range_write
    module subroutine range_int_write (object, unit)
      class(range_int_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine range_int_write
    module subroutine range_real_write (object, unit)
      class(range_real_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine range_real_write
    module subroutine range_init (range, pn)
      class(range_t), intent(out) :: range
      type(parse_node_t), intent(in), target :: pn
    end subroutine range_init
    module subroutine range_create_value_node (range)
      class(range_t), intent(inout) :: range
    end subroutine range_create_value_node
    module subroutine range_compile (range, global)
      class(range_t), intent(inout) :: range
      type(rt_data_t), intent(in), target :: global
    end subroutine range_compile
    module subroutine range_int_evaluate (range)
      class(range_int_t), intent(inout) :: range
    end subroutine range_int_evaluate
    module subroutine range_real_evaluate (range)
      class(range_real_t), intent(inout) :: range
    end subroutine range_real_evaluate
    module function range_get_n_iterations (range) result (n)
      class(range_t), intent(in) :: range
      integer :: n
    end function range_get_n_iterations
    module subroutine range_int_set_value (range, i)
      class(range_int_t), intent(inout) :: range
      integer, intent(in) :: i
    end subroutine range_int_set_value
    module subroutine range_real_set_value (range, i)
      class(range_real_t), intent(inout) :: range
      integer, intent(in) :: i
    end subroutine range_real_set_value
    recursive module subroutine cmd_scan_final (cmd)
      class(cmd_scan_t), intent(inout) :: cmd
    end subroutine cmd_scan_final
    module subroutine cmd_scan_write (cmd, unit, indent)
      class(cmd_scan_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_scan_write
    recursive module subroutine cmd_scan_execute (cmd, global)
      class(cmd_scan_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_scan_execute
    recursive module subroutine cmd_if_final (cmd)
      class(cmd_if_t), intent(inout) :: cmd
    end subroutine cmd_if_final
    module subroutine cmd_if_write (cmd, unit, indent)
      class(cmd_if_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_if_write
    recursive module subroutine cmd_if_compile (cmd, global)
      class(cmd_if_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_if_compile
    recursive module subroutine cmd_if_execute (cmd, global)
      class(cmd_if_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_if_execute
    module subroutine cmd_include_final (cmd)
      class(cmd_include_t), intent(inout) :: cmd
    end subroutine cmd_include_final
    module subroutine cmd_include_write (cmd, unit, indent)
      class(cmd_include_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_include_write
    module subroutine cmd_include_compile (cmd, global)
      class(cmd_include_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_include_compile
    module subroutine cmd_include_execute (cmd, global)
      class(cmd_include_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_include_execute
    module subroutine cmd_export_write (cmd, unit, indent)
      class(cmd_export_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_export_write
    module subroutine cmd_export_compile (cmd, global)
      class(cmd_export_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_export_compile
    module subroutine cmd_export_execute (cmd, global)
      class(cmd_export_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_export_execute
    module subroutine cmd_quit_write (cmd, unit, indent)
      class(cmd_quit_t), intent(in) :: cmd
      integer, intent(in), optional :: unit, indent
    end subroutine cmd_quit_write
    module subroutine cmd_quit_compile (cmd, global)
      class(cmd_quit_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_quit_compile
    module subroutine cmd_quit_execute (cmd, global)
      class(cmd_quit_t), intent(inout) :: cmd
      type(rt_data_t), intent(inout), target :: global
    end subroutine cmd_quit_execute
    recursive module subroutine command_list_write (cmd_list, unit, indent)
      class(command_list_t), intent(in) :: cmd_list
      integer, intent(in), optional :: unit, indent
    end subroutine command_list_write
    module subroutine command_list_append (cmd_list, command)
      class(command_list_t), intent(inout) :: cmd_list
      class(command_t), intent(inout), pointer :: command
    end subroutine command_list_append
    recursive module subroutine command_list_final (cmd_list)
      class(command_list_t), intent(inout) :: cmd_list
    end subroutine command_list_final
    recursive module subroutine command_list_execute (cmd_list, global)
      class(command_list_t), intent(in) :: cmd_list
      type(rt_data_t), intent(inout), target :: global
    end subroutine command_list_execute
    module subroutine syntax_cmd_list_init ()
    end subroutine syntax_cmd_list_init
    module subroutine syntax_cmd_list_final ()
    end subroutine syntax_cmd_list_final
    module subroutine syntax_cmd_list_write (unit)
      integer, intent(in), optional :: unit
    end subroutine syntax_cmd_list_write
    module subroutine lexer_init_cmd_list (lexer, parent_lexer)
      type(lexer_t), intent(out) :: lexer
      type(lexer_t), intent(in), optional, target :: parent_lexer
    end subroutine lexer_init_cmd_list
  end interface

contains

  subroutine dispatch_command (command, pn)
    class(command_t), intent(inout), pointer :: command
    type(parse_node_t), intent(in), target :: pn
    select case (char (parse_node_get_rule_key (pn)))
    case ("cmd_model")
       allocate (cmd_model_t :: command)
    case ("cmd_library")
       allocate (cmd_library_t :: command)
    case ("cmd_process")
       allocate (cmd_process_t :: command)
    case ("cmd_nlo")
       allocate (cmd_nlo_t :: command)
    case ("cmd_compile")
       allocate (cmd_compile_t :: command)
    case ("cmd_exec")
       allocate (cmd_exec_t :: command)
     case ("cmd_num", "cmd_complex", "cmd_real", "cmd_int", &
           "cmd_log_decl", "cmd_log", "cmd_string", "cmd_string_decl", &
           "cmd_alias", "cmd_result")
       allocate (cmd_var_t :: command)
    case ("cmd_slha")
       allocate (cmd_slha_t :: command)
    case ("cmd_show")
       allocate (cmd_show_t :: command)
    case ("cmd_clear")
       allocate (cmd_clear_t :: command)
    case ("cmd_expect")
       allocate (cmd_expect_t :: command)
    case ("cmd_beams")
       allocate (cmd_beams_t :: command)
    case ("cmd_beams_pol_density")
       allocate (cmd_beams_pol_density_t :: command)
    case ("cmd_beams_pol_fraction")
       allocate (cmd_beams_pol_fraction_t :: command)
    case ("cmd_beams_momentum")
       allocate (cmd_beams_momentum_t :: command)
    case ("cmd_beams_theta")
       allocate (cmd_beams_theta_t :: command)
    case ("cmd_beams_phi")
       allocate (cmd_beams_phi_t :: command)
    case ("cmd_cuts")
       allocate (cmd_cuts_t :: command)
    case ("cmd_scale")
       allocate (cmd_scale_t :: command)
    case ("cmd_fac_scale")
       allocate (cmd_fac_scale_t :: command)
    case ("cmd_ren_scale")
       allocate (cmd_ren_scale_t :: command)
    case ("cmd_weight")
       allocate (cmd_weight_t :: command)
    case ("cmd_selection")
       allocate (cmd_selection_t :: command)
    case ("cmd_reweight")
       allocate (cmd_reweight_t :: command)
    case ("cmd_iterations")
       allocate (cmd_iterations_t :: command)
    case ("cmd_integrate")
       allocate (cmd_integrate_t :: command)
    case ("cmd_observable")
       allocate (cmd_observable_t :: command)
    case ("cmd_histogram")
       allocate (cmd_histogram_t :: command)
    case ("cmd_plot")
       allocate (cmd_plot_t :: command)
    case ("cmd_graph")
       allocate (cmd_graph_t :: command)
    case ("cmd_record")
       allocate (cmd_record_t :: command)
    case ("cmd_analysis")
       allocate (cmd_analysis_t :: command)
    case ("cmd_alt_setup")
       allocate (cmd_alt_setup_t :: command)
    case ("cmd_unstable")
       allocate (cmd_unstable_t :: command)
    case ("cmd_stable")
       allocate (cmd_stable_t :: command)
    case ("cmd_polarized")
       allocate (cmd_polarized_t :: command)
    case ("cmd_unpolarized")
       allocate (cmd_unpolarized_t :: command)
    case ("cmd_sample_format")
       allocate (cmd_sample_format_t :: command)
    case ("cmd_simulate")
       allocate (cmd_simulate_t :: command)
    case ("cmd_rescan")
       allocate (cmd_rescan_t :: command)
    case ("cmd_write_analysis")
       allocate (cmd_write_analysis_t :: command)
    case ("cmd_compile_analysis")
       allocate (cmd_compile_analysis_t :: command)
    case ("cmd_open_out")
       allocate (cmd_open_out_t :: command)
    case ("cmd_close_out")
       allocate (cmd_close_out_t :: command)
    case ("cmd_printf")
       allocate (cmd_printf_t :: command)
    case ("cmd_scan")
       allocate (cmd_scan_t :: command)
    case ("cmd_if")
       allocate (cmd_if_t :: command)
    case ("cmd_include")
       allocate (cmd_include_t :: command)
    case ("cmd_export")
       allocate (cmd_export_t :: command)
    case ("cmd_quit")
       allocate (cmd_quit_t :: command)
    case default
       print *, char (parse_node_get_rule_key (pn))
       call msg_bug ("Command not implemented")
    end select
    command%pn => pn
  end subroutine dispatch_command

  recursive subroutine cmd_scan_compile (cmd, global)
    class(cmd_scan_t), intent(inout) :: cmd
    type(rt_data_t), intent(inout), target :: global
    type(var_list_t), pointer :: var_list
    type(parse_node_t), pointer :: pn_var, pn_body, pn_body_first
    type(parse_node_t), pointer :: pn_decl, pn_name
    type(parse_node_t), pointer :: pn_arg, pn_scan_cmd, pn_rhs
    type(parse_node_t), pointer :: pn_decl_single, pn_var_single
    type(syntax_rule_t), pointer :: var_rule_decl, var_rule
    type(string_t) :: key
    integer :: var_type
    integer :: i
    if (debug_on) call msg_debug (D_CORE, "cmd_scan_compile")
    if (debug_active (D_CORE))  call parse_node_write_rec (cmd%pn)
    pn_var => parse_node_get_sub_ptr (cmd%pn, 2)
    pn_body => parse_node_get_next_ptr (pn_var)
    if (associated (pn_body)) then
       pn_body_first => parse_node_get_sub_ptr (pn_body)
    else
       pn_body_first => null ()
    end if
    key = parse_node_get_rule_key (pn_var)
    select case (char (key))
    case ("scan_num")
       pn_name => parse_node_get_sub_ptr (pn_var)
       cmd%name = parse_node_get_string (pn_name)
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, var_str ("cmd_num"))
       pn_arg => parse_node_get_next_ptr (pn_name, 2)
    case ("scan_int")
       pn_name => parse_node_get_sub_ptr (pn_var, 2)
       cmd%name = parse_node_get_string (pn_name)
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, var_str ("cmd_int"))
       pn_arg => parse_node_get_next_ptr (pn_name, 2)
    case ("scan_real")
       pn_name => parse_node_get_sub_ptr (pn_var, 2)
       cmd%name = parse_node_get_string (pn_name)
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, var_str ("cmd_real"))
       pn_arg => parse_node_get_next_ptr (pn_name, 2)
    case ("scan_complex")
       pn_name => parse_node_get_sub_ptr (pn_var, 2)
       cmd%name = parse_node_get_string (pn_name)
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, var_str("cmd_complex"))
       pn_arg => parse_node_get_next_ptr (pn_name, 2)
    case ("scan_alias")
       pn_name => parse_node_get_sub_ptr (pn_var, 2)
       cmd%name = parse_node_get_string (pn_name)
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, var_str ("cmd_alias"))
       pn_arg => parse_node_get_next_ptr (pn_name, 2)
    case ("scan_string_decl")
       pn_decl => parse_node_get_sub_ptr (pn_var, 2)
       pn_name => parse_node_get_sub_ptr (pn_decl, 2)
       cmd%name = parse_node_get_string (pn_name)
       var_rule_decl => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_string"))
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_string_decl"))
       pn_arg => parse_node_get_next_ptr (pn_name, 2)
    case ("scan_log_decl")
       pn_decl => parse_node_get_sub_ptr (pn_var, 2)
       pn_name => parse_node_get_sub_ptr (pn_decl, 2)
       cmd%name = parse_node_get_string (pn_name)
       var_rule_decl => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_log"))
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_log_decl"))
       pn_arg => parse_node_get_next_ptr (pn_name, 2)
    case ("scan_cuts")
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_cuts"))
       cmd%name = "cuts"
       pn_arg => parse_node_get_sub_ptr (pn_var, 3)
    case ("scan_weight")
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_weight"))
       cmd%name = "weight"
       pn_arg => parse_node_get_sub_ptr (pn_var, 3)
    case ("scan_scale")
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_scale"))
       cmd%name = "scale"
       pn_arg => parse_node_get_sub_ptr (pn_var, 3)
    case ("scan_ren_scale")
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_ren_scale"))
       cmd%name = "renormalization_scale"
       pn_arg => parse_node_get_sub_ptr (pn_var, 3)
    case ("scan_fac_scale")
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_fac_scale"))
       cmd%name = "factorization_scale"
       pn_arg => parse_node_get_sub_ptr (pn_var, 3)
    case ("scan_selection")
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_selection"))
       cmd%name = "selection"
       pn_arg => parse_node_get_sub_ptr (pn_var, 3)
    case ("scan_reweight")
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_reweight"))
       cmd%name = "reweight"
       pn_arg => parse_node_get_sub_ptr (pn_var, 3)
    case ("scan_analysis")
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_analysis"))
       cmd%name = "analysis"
       pn_arg => parse_node_get_sub_ptr (pn_var, 3)
    case ("scan_model")
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_model"))
       cmd%name = "model"
       pn_arg => parse_node_get_sub_ptr (pn_var, 3)
    case ("scan_library")
       var_rule => syntax_get_rule_ptr (syntax_cmd_list, &
            var_str ("cmd_library"))
       cmd%name = "library"
       pn_arg => parse_node_get_sub_ptr (pn_var, 3)
    case default
       call msg_bug ("scan: case '" // char (key) // "' not implemented")
    end select
    if (associated (pn_arg)) then
       cmd%n_values = parse_node_get_n_sub (pn_arg)
    end if
    var_list => global%get_var_list_ptr ()
    allocate (cmd%scan_cmd (cmd%n_values))
    select case (char (key))
    case ("scan_num")
       var_type = &
            var_list%get_type (cmd%name)
       select case (var_type)
       case (V_INT)
          allocate (range_int_t :: cmd%range (cmd%n_values))
       case (V_REAL)
          allocate (range_real_t :: cmd%range (cmd%n_values))
       case (V_CMPLX)
          call msg_fatal ("scan over complex variable not implemented")
       case (V_NONE)
          call msg_fatal ("scan: variable '" // char (cmd%name) //"' undefined")
       case default
          call msg_bug ("scan: impossible variable type")
       end select
    case ("scan_int")
       allocate (range_int_t :: cmd%range (cmd%n_values))
    case ("scan_real")
       allocate (range_real_t :: cmd%range (cmd%n_values))
    case ("scan_complex")
       call msg_fatal ("scan over complex variable not implemented")
    end select
    i = 1
    if (associated (pn_arg)) then
       pn_rhs => parse_node_get_sub_ptr (pn_arg)
    else
       pn_rhs => null ()
    end if
    do while (associated (pn_rhs))
       allocate (pn_scan_cmd)
       call parse_node_create_branch (pn_scan_cmd, &
            syntax_get_rule_ptr (syntax_cmd_list, var_str ("command_list")))
       allocate (pn_var_single)
       pn_var_single = pn_var
       call parse_node_replace_rule (pn_var_single, var_rule)
       select case (char (key))
       case ("scan_num", "scan_int", "scan_real", &
            "scan_complex", "scan_alias", &
            "scan_cuts", "scan_weight", &
            "scan_scale", "scan_ren_scale", "scan_fac_scale", &
            "scan_selection", "scan_reweight", "scan_analysis", &
            "scan_model", "scan_library")
          if (allocated (cmd%range)) then
             call cmd%range(i)%init (pn_rhs)
             call parse_node_replace_last_sub &
                  (pn_var_single, cmd%range(i)%pn_expr)
          else
             call parse_node_replace_last_sub (pn_var_single, pn_rhs)
          end if
       case ("scan_string_decl", "scan_log_decl")
          allocate (pn_decl_single)
          pn_decl_single = pn_decl
          call parse_node_replace_rule (pn_decl_single, var_rule_decl)
          call parse_node_replace_last_sub (pn_decl_single, pn_rhs)
          call parse_node_freeze_branch (pn_decl_single)
          call parse_node_replace_last_sub (pn_var_single, pn_decl_single)
       case default
          call msg_bug ("scan: case '" // char (key)  &
               // "' broken")
       end select
       call parse_node_freeze_branch (pn_var_single)
       call parse_node_append_sub (pn_scan_cmd, pn_var_single)
       call parse_node_append_sub (pn_scan_cmd, pn_body_first)
       call parse_node_freeze_branch (pn_scan_cmd)
       cmd%scan_cmd(i)%ptr => pn_scan_cmd
       i = i + 1
       pn_rhs => parse_node_get_next_ptr (pn_rhs)
    end do
    if (debug_active (D_CORE)) then
       do i = 1, cmd%n_values
          print *, "scan command ", i
          call parse_node_write_rec (cmd%scan_cmd(i)%ptr)
          if (allocated (cmd%range))  call cmd%range(i)%write ()
       end do
       print *, "original"
       call parse_node_write_rec (cmd%pn)
    end if
  end subroutine cmd_scan_compile

  recursive subroutine command_list_compile (cmd_list, pn, global)
    class(command_list_t), intent(inout), target :: cmd_list
    type(parse_node_t), intent(in), target :: pn
    type(rt_data_t), intent(inout), target :: global
    type(parse_node_t), pointer :: pn_cmd
    class(command_t), pointer :: command
    integer :: i
    pn_cmd => parse_node_get_sub_ptr (pn)
    do i = 1, parse_node_get_n_sub (pn)
       call dispatch_command (command, pn_cmd)
       call command%compile (global)
       call cmd_list%append (command)
       call terminate_now_if_signal ()
       pn_cmd => parse_node_get_next_ptr (pn_cmd)
    end do
  end subroutine command_list_compile


end module commands
