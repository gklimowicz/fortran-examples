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

module processes_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use format_utils, only: write_separator
  use constants, only: TWOPI4
  use physics_defs, only: CONV
  use os_interface
  use sm_qcd
  use lorentz
  use pdg_arrays
  use model_data
  use models
  use var_base, only: vars_t
  use variables, only: var_list_t
  use model_testbed, only: prepare_model
  use particle_specifiers, only: new_prt_spec
  use flavors
  use interactions, only: reset_interaction_counter
  use particles
  use rng_base
  use mci_base
  use mci_none, only: mci_none_t
  use mci_midpoint
  use sf_mappings
  use sf_base
  use phs_base
  use phs_single
  use phs_forests, only: syntax_phs_forest_init, syntax_phs_forest_final
  use phs_wood, only: phs_wood_config_t
  use resonances, only: resonance_history_set_t
  use process_constants
  use prc_core_def, only: prc_core_def_t
  use prc_core
  use prc_test, only: prc_test_create_library
  use prc_template_me, only: template_me_def_t
  use process_libraries
  use prc_test_core
  use pdf, only: pdf_data_t

  use process_counter
  use process_config, only: process_term_t
  use process, only: process_t
  use instances, only: process_instance_t, process_instance_hook_t

  use rng_base_ut, only: rng_test_factory_t
  use sf_base_ut, only: sf_test_data_t
  use mci_base_ut, only: mci_test_t
  use phs_base_ut, only: phs_test_config_t

  implicit none
  private

  public :: prepare_test_process
  public :: cleanup_test_process

  public :: processes_1
  public :: processes_2
  public :: processes_3
  public :: processes_4
  public :: processes_7
  public :: processes_8
  public :: processes_9
  public :: processes_10
  public :: processes_11
  public :: processes_12
  public :: processes_13
  public :: processes_14
  public :: processes_15
  public :: processes_16
  public :: processes_17
  public :: processes_18
  public :: processes_19

  type, extends(process_instance_hook_t) :: process_instance_hook_test_t
    integer :: unit
    character(len=15) :: name
  contains
    procedure :: init => process_instance_hook_test_init
    procedure :: final => process_instance_hook_test_final
    procedure :: evaluate => process_instance_hook_test_evaluate
  end type process_instance_hook_test_t


contains

  subroutine processes_1 (u)
    integer, intent(in) :: u
    type(process_t) :: process

    write (u, "(A)")  "* Test output: processes_1"
    write (u, "(A)")  "*   Purpose: display an empty process object"
    write (u, "(A)")

    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_1"

  end subroutine processes_1

  subroutine processes_2 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable :: process
    class(mci_t), allocatable :: mci_template
    class(phs_config_t), allocatable :: phs_config_template

    write (u, "(A)")  "* Test output: processes_2"
    write (u, "(A)")  "*   Purpose: initialize a simple process object"
    write (u, "(A)")

    write (u, "(A)")  "* Build and load a test library with one process"
    write (u, "(A)")

    libname = "processes2"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib)

    write (u, "(A)")  "* Initialize a process object"
    write (u, "(A)")

    call model%init_test ()

    allocate (process)
    call process%init (procname, lib, os_data, model)
    call process%set_run_id (var_str ("run_2"))
    call process%setup_test_cores ()

    allocate (phs_test_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    call process%setup_mci (dispatch_mci_empty)

    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_2"

  end subroutine processes_2

  subroutine processes_3 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable :: process
    class(phs_config_t), allocatable :: phs_config_template
    type(process_constants_t) :: data
    type(vector4_t), dimension(:), allocatable :: p

    write (u, "(A)")  "* Test output: processes_3"
    write (u, "(A)")  "*   Purpose: create a process &
         &and compute a matrix element"
    write (u, "(A)")

    write (u, "(A)")  "* Build and load a test library with one process"
    write (u, "(A)")

    libname = "processes3"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib)

    call model%init_test ()

    allocate (process)
    call process%init (procname, lib, os_data, model)
    call process%setup_test_cores ()

    allocate (phs_test_config_t :: phs_config_template)
    call process%init_components (phs_config_template)
    call process%setup_mci (dispatch_mci_test3)

    write (u, "(A)")  "* Return the number of process components"
    write (u, "(A)")

    write (u, "(A,I0)")  "n_components = ", process%get_n_components ()

    write (u, "(A)")
    write (u, "(A)")  "* Return the number of flavor states"
    write (u, "(A)")

    data = process%get_constants (1)

    write (u, "(A,I0)")  "n_flv(1) = ", data%n_flv

    write (u, "(A)")
    write (u, "(A)")  "* Return the first flavor state"
    write (u, "(A)")

    write (u, "(A,4(1x,I0))")  "flv_state(1) =", data%flv_state (:,1)

    write (u, "(A)")
    write (u, "(A)")  "* Set up kinematics &
         &[arbitrary, the matrix element is constant]"

    allocate (p (4))

    write (u, "(A)")
    write (u, "(A)")  "* Retrieve the matrix element"
    write (u, "(A)")


    write (u, "(A,F5.3,' + ',F5.3,' I')")  "me (1, p, 1, 1, 1) = ", &
         process%compute_amplitude (1, 1, 1, p, 1, 1, 1)


    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_3"

  end subroutine processes_3

  subroutine processes_4 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable, target :: process
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts
    type(process_instance_t), allocatable, target :: process_instance
    type(particle_set_t) :: pset

    write (u, "(A)")  "* Test output: processes_4"
    write (u, "(A)")  "*   Purpose: create a process &
         &and fill a process instance"
    write (u, "(A)")

    write (u, "(A)")  "* Build and initialize a test process"
    write (u, "(A)")

    libname = "processes4"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib)

    call reset_interaction_counter ()

    call model%init_test ()

    allocate (process)
    call process%init (procname, lib, os_data, model)

    call process%setup_test_cores ()
    allocate (phs_test_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    write (u, "(A)")  "* Prepare a trivial beam setup"
    write (u, "(A)")

    sqrts = 1000
    call process%setup_beams_sqrts (sqrts, i_core = 1)
    call process%configure_phs ()
    call process%setup_mci (dispatch_mci_empty)

    write (u, "(A)")  "* Complete process initialization"
    write (u, "(A)")

    call process%setup_terms ()
    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Create a process instance"
    write (u, "(A)")

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Inject a set of random numbers"
    write (u, "(A)")

    call process_instance%choose_mci (1)
    call process_instance%set_mcpar ([0._default, 0._default])
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Set up hard kinematics"
    write (u, "(A)")

    call process_instance%select_channel (1)
    call process_instance%compute_seed_kinematics ()
    call process_instance%compute_hard_kinematics ()
    call process_instance%compute_eff_kinematics ()
    call process_instance%evaluate_expressions ()
    call process_instance%compute_other_channels ()

    write (u, "(A)")  "* Evaluate matrix element and square"
    write (u, "(A)")

    call process_instance%evaluate_trace ()
    call process_instance%write (u)

    call process_instance%get_trace (pset, 1)
    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Particle content:"
    write (u, "(A)")

    call write_separator (u)
    call pset%write (u)
    call write_separator (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover process instance"
    write (u, "(A)")

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%choose_mci (1)
    call process_instance%set_trace (pset, 1, check_match = .false.)

    call process_instance%activate ()
    process_instance%evaluation_status = STAT_EFF_KINEMATICS
    call process_instance%recover_hard_kinematics (i_term = 1)
    call process_instance%recover_seed_kinematics (i_term = 1)
    call process_instance%select_channel (1)
    call process_instance%recover_mcpar (i_term = 1)

    call process_instance%compute_seed_kinematics (skip_term = 1)
    call process_instance%compute_hard_kinematics (skip_term = 1)
    call process_instance%compute_eff_kinematics (skip_term = 1)

    call process_instance%evaluate_expressions ()
    call process_instance%compute_other_channels (skip_term = 1)
    call process_instance%evaluate_trace ()
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call pset%final ()
    call process_instance%final ()
    deallocate (process_instance)

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_4"

  end subroutine processes_4

  subroutine processes_7 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable, target :: process
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    type(sf_config_t), dimension(:), allocatable :: sf_config
    type(sf_channel_t), dimension(2) :: sf_channel

    write (u, "(A)")  "* Test output: processes_7"
    write (u, "(A)")  "*   Purpose: initialize a process with &
         &structure functions"
    write (u, "(A)")

    write (u, "(A)")  "* Build and initialize a process object"
    write (u, "(A)")

    libname = "processes7"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib)

    call model%init_test ()

    allocate (process)
    call process%init (procname, lib, os_data, model)

    call process%setup_test_cores ()
    allocate (phs_test_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    write (u, "(A)")  "* Set beam, structure functions, and mappings"
    write (u, "(A)")

    sqrts = 1000
    call process%setup_beams_sqrts (sqrts, i_core = 1)
    call process%configure_phs ()

    pdg_in = 25
    allocate (sf_test_data_t :: data)
    select type (data)
    type is (sf_test_data_t)
       call data%init (process%get_model_ptr (), pdg_in)
    end select

    allocate (sf_config (2))
    call sf_config(1)%init ([1], data)
    call sf_config(2)%init ([2], data)
    call process%init_sf_chain (sf_config)
    deallocate (sf_config)

    call process%test_allocate_sf_channels (3)

    call sf_channel(1)%init (2)
    call sf_channel(1)%activate_mapping ([1,2])
    call process%set_sf_channel (2, sf_channel(1))

    call sf_channel(2)%init (2)
    call sf_channel(2)%set_s_mapping ([1,2])
    call process%set_sf_channel (3, sf_channel(2))

    call process%setup_mci (dispatch_mci_empty)

    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_7"

  end subroutine processes_7

  subroutine processes_8 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable, target :: process
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts
    type(process_instance_t), allocatable, target :: process_instance
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    type(sf_config_t), dimension(:), allocatable :: sf_config
    type(sf_channel_t) :: sf_channel
    type(particle_set_t) :: pset

    write (u, "(A)")  "* Test output: processes_8"
    write (u, "(A)")  "*   Purpose: evaluate a process with &
         &structure functions"
    write (u, "(A)")

    write (u, "(A)")  "* Build and initialize a process object"
    write (u, "(A)")

    libname = "processes8"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib)

    call reset_interaction_counter ()

    call model%init_test ()

    allocate (process)
    call process%init (procname, lib, os_data, model)

    call process%setup_test_cores ()
    allocate (phs_test_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    write (u, "(A)")  "* Set beam, structure functions, and mappings"
    write (u, "(A)")

    sqrts = 1000
    call process%setup_beams_sqrts (sqrts, i_core = 1)

    pdg_in = 25
    allocate (sf_test_data_t :: data)
    select type (data)
    type is (sf_test_data_t)
       call data%init (process%get_model_ptr (), pdg_in)
    end select

    allocate (sf_config (2))
    call sf_config(1)%init ([1], data)
    call sf_config(2)%init ([2], data)
    call process%init_sf_chain (sf_config)
    deallocate (sf_config)

    call process%configure_phs ()

    call process%test_allocate_sf_channels (1)

    call sf_channel%init (2)
    call sf_channel%activate_mapping ([1,2])
    call process%set_sf_channel (1, sf_channel)

    write (u, "(A)")  "* Complete process initialization"
    write (u, "(A)")

    call process%setup_mci (dispatch_mci_empty)
    call process%setup_terms ()

    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Create a process instance"
    write (u, "(A)")

    allocate (process_instance)
    call process_instance%init (process)

    write (u, "(A)")  "* Set up kinematics and evaluate"
    write (u, "(A)")

    call process_instance%choose_mci (1)
    call process_instance%evaluate_sqme (1, &
         [0.8_default, 0.8_default, 0.1_default, 0.2_default])
    call process_instance%write (u)

    call process_instance%get_trace (pset, 1)
    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Particle content:"
    write (u, "(A)")

    call write_separator (u)
    call pacify (pset)
    call pset%write (u)
    call write_separator (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover process instance"
    write (u, "(A)")

    call reset_interaction_counter (2)

    allocate (process_instance)
    call process_instance%init (process)

    call process_instance%choose_mci (1)
    call process_instance%set_trace (pset, 1, check_match = .false.)
    call process_instance%recover &
         (channel = 1, i_term = 1, update_sqme = .true., recover_phs = .true.)
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call pset%final ()

    call process_instance%final ()
    deallocate (process_instance)

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_8"

  end subroutine processes_8

  subroutine processes_9 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable, target :: process
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts
    type(process_instance_t), allocatable, target :: process_instance
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    type(sf_config_t), dimension(:), allocatable :: sf_config
    type(sf_channel_t) :: sf_channel
    real(default), dimension(4) :: x_saved
    type(particle_set_t) :: pset

    write (u, "(A)")  "* Test output: processes_9"
    write (u, "(A)")  "*   Purpose: evaluate a process with &
         &structure functions"
    write (u, "(A)")  "*            in a multi-channel configuration"
    write (u, "(A)")

    write (u, "(A)")  "* Build and initialize a process object"
    write (u, "(A)")

    libname = "processes9"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib)

    call reset_interaction_counter ()

    call model%init_test ()

    allocate (process)
    call process%init (procname, lib, os_data, model)

    call process%setup_test_cores ()
    allocate (phs_test_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    write (u, "(A)")  "* Set beam, structure functions, and mappings"
    write (u, "(A)")

    sqrts = 1000
    call process%setup_beams_sqrts (sqrts, i_core = 1)

    pdg_in = 25
    allocate (sf_test_data_t :: data)
    select type (data)
    type is (sf_test_data_t)
       call data%init (process%get_model_ptr (), pdg_in)
    end select

    allocate (sf_config (2))
    call sf_config(1)%init ([1], data)
    call sf_config(2)%init ([2], data)
    call process%init_sf_chain (sf_config)
    deallocate (sf_config)

    call process%configure_phs ()

    call process%test_allocate_sf_channels (2)

    call sf_channel%init (2)
    call process%set_sf_channel (1, sf_channel)

    call sf_channel%init (2)
    call sf_channel%activate_mapping ([1,2])
    call process%set_sf_channel (2, sf_channel)

    call process%test_set_component_sf_channel ([1, 2])

    write (u, "(A)")  "* Complete process initialization"
    write (u, "(A)")

    call process%setup_mci (dispatch_mci_empty)
    call process%setup_terms ()

    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Create a process instance"
    write (u, "(A)")

    allocate (process_instance)
    call process_instance%init (process)

    write (u, "(A)")  "* Set up kinematics in channel 1 and evaluate"
    write (u, "(A)")

    call process_instance%choose_mci (1)
    call process_instance%evaluate_sqme (1, &
         [0.8_default, 0.8_default, 0.1_default, 0.2_default])
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Extract MC input parameters"
    write (u, "(A)")

    write (u, "(A)")  "Channel 1:"
    call process_instance%get_mcpar (1, x_saved)
    write (u, "(2x,9(1x,F7.5))")  x_saved

    write (u, "(A)")  "Channel 2:"
    call process_instance%get_mcpar (2, x_saved)
    write (u, "(2x,9(1x,F7.5))")  x_saved

    write (u, "(A)")
    write (u, "(A)")  "* Set up kinematics in channel 2 and evaluate"
    write (u, "(A)")

    call process_instance%evaluate_sqme (2, x_saved)
    call process_instance%write (u)

    call process_instance%get_trace (pset, 1)
    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Recover process instance for channel 2"
    write (u, "(A)")

    call reset_interaction_counter (2)

    allocate (process_instance)
    call process_instance%init (process)

    call process_instance%choose_mci (1)
    call process_instance%set_trace (pset, 1, check_match = .false.)
    call process_instance%recover &
         (channel = 2, i_term = 1, update_sqme = .true., recover_phs = .true.)
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call pset%final ()

    call process_instance%final ()
    deallocate (process_instance)

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_9"

  end subroutine processes_9

  subroutine processes_10 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable, target :: process
    class(mci_t), pointer :: mci
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts
    type(process_instance_t), allocatable, target :: process_instance

    write (u, "(A)")  "* Test output: processes_10"
    write (u, "(A)")  "*   Purpose: generate events for a process without &
         &structure functions"
    write (u, "(A)")  "*            in a multi-channel configuration"
    write (u, "(A)")

    write (u, "(A)")  "* Build and initialize a process object"
    write (u, "(A)")

    libname = "processes10"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib)

    call reset_interaction_counter ()

    call model%init_test ()

    allocate (process)
    call process%init (procname, lib, os_data, model)

    call process%setup_test_cores ()
    allocate (phs_test_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    write (u, "(A)")  "* Prepare a trivial beam setup"
    write (u, "(A)")

    sqrts = 1000
    call process%setup_beams_sqrts (sqrts, i_core = 1)
    call process%configure_phs ()

    call process%setup_mci (dispatch_mci_test10)

    write (u, "(A)")  "* Complete process initialization"
    write (u, "(A)")

    call process%setup_terms ()
    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Create a process instance"
    write (u, "(A)")

    allocate (process_instance)
    call process_instance%init (process)

    write (u, "(A)")  "* Generate weighted event"
    write (u, "(A)")

    call process%test_get_mci_ptr (mci)
    select type (mci)
    type is (mci_test_t)
       ! This ensures that the next 'random' numbers are 0.3, 0.5, 0.7
       call mci%rng%init (3)
       ! Include the constant PHS factor in the stored maximum of the integrand
       call mci%set_max_factor (conv * twopi4 &
            / (2 * sqrt (lambda (sqrts **2, 125._default**2, 125._default**2))))
    end select

    call process_instance%generate_weighted_event (1)
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate unweighted event"
    write (u, "(A)")

    call process_instance%generate_unweighted_event (1)
    call process%test_get_mci_ptr (mci)
    select type (mci)
    type is (mci_test_t)
       write (u, "(A,I0)")  " Success in try ", mci%tries
       write (u, "(A)")
    end select

    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call process_instance%final ()
    deallocate (process_instance)

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_10"

  end subroutine processes_10

  subroutine processes_11 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable, target :: process
    class(mci_t), allocatable :: mci_template
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts
    type(process_instance_t), allocatable, target :: process_instance

    write (u, "(A)")  "* Test output: processes_11"
    write (u, "(A)")  "*   Purpose: integrate a process without &
         &structure functions"
    write (u, "(A)")  "*            in a multi-channel configuration"
    write (u, "(A)")

    write (u, "(A)")  "* Build and initialize a process object"
    write (u, "(A)")

    libname = "processes11"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib)

    call reset_interaction_counter ()

    call model%init_test ()

    allocate (process)
    call process%init (procname, lib, os_data, model)

    call process%setup_test_cores ()

    allocate (phs_test_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    write (u, "(A)")  "* Prepare a trivial beam setup"
    write (u, "(A)")

    sqrts = 1000
    call process%setup_beams_sqrts (sqrts, i_core = 1)
    call process%configure_phs ()

    call process%setup_mci (dispatch_mci_test10)

    write (u, "(A)")  "* Complete process initialization"
    write (u, "(A)")

    call process%setup_terms ()
    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Create a process instance"
    write (u, "(A)")

    allocate (process_instance)
    call process_instance%init (process)

    write (u, "(A)")  "* Integrate with default test parameters"
    write (u, "(A)")

    call process_instance%integrate (1, n_it=1, n_calls=10000)
    call process%final_integration (1)

    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A,ES13.7)")  " Integral divided by phs factor = ", &
         process%get_integral (1) &
         / process_instance%kin(1)%phs_factor

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call process_instance%final ()
    deallocate (process_instance)

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_11"

  end subroutine processes_11

  subroutine processes_12 (u)
    integer, intent(in) :: u
    type(process_t), allocatable, target :: process
    type(process_instance_t), allocatable, target :: process_instance
    type(particle_set_t) :: pset
    type(model_data_t), target :: model

    write (u, "(A)")  "* Test output: processes_12"
    write (u, "(A)")  "*   Purpose: generate a complete partonic event"
    write (u, "(A)")

    call model%init_test ()

    write (u, "(A)")  "* Build and initialize process and process instance &
         &and generate event"
    write (u, "(A)")

    allocate (process)
    allocate (process_instance)
    call prepare_test_process (process, process_instance, model, &
         run_id = var_str ("run_12"))
    call process_instance%setup_event_data (i_core = 1)

    call process%prepare_simulation (1)
    call process_instance%init_simulation (1)
    call process_instance%generate_weighted_event (1)
    call process_instance%evaluate_event_data ()

    call process_instance%write (u)

    call process_instance%get_trace (pset, 1)

    call process_instance%final_simulation (1)
    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Recover kinematics and recalculate"
    write (u, "(A)")

    call reset_interaction_counter (2)

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data ()

    call process_instance%choose_mci (1)
    call process_instance%set_trace (pset, 1, check_match = .false.)
    call process_instance%recover &
         (channel = 1, i_term = 1, update_sqme = .true., recover_phs = .true.)

    call process_instance%recover_event ()
    call process_instance%evaluate_event_data ()

    call process_instance%write (u)


    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call cleanup_test_process (process, process_instance)
    deallocate (process_instance)
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_12"

  end subroutine processes_12

  subroutine processes_13 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(model_data_t), target :: model
    type(process_term_t) :: term
    class(prc_core_t), allocatable :: core

    write (u, "(A)")  "* Test output: processes_13"
    write (u, "(A)")  "*   Purpose: initialized a colored interaction"
    write (u, "(A)")

    write (u, "(A)")  "* Set up a process constants block"
    write (u, "(A)")

    call os_data%init ()
    call model%init_sm_test ()

    allocate (test_t :: core)

    associate (data => term%data)
      data%n_in = 2
      data%n_out = 3
      data%n_flv = 2
      data%n_hel = 2
      data%n_col = 2
      data%n_cin = 2

      allocate (data%flv_state (5, 2))
      data%flv_state (:,1) = [ 1, 21, 1, 21, 21]
      data%flv_state (:,2) = [ 2, 21, 2, 21, 21]

      allocate (data%hel_state (5, 2))
      data%hel_state (:,1) = [1, 1, 1, 1, 0]
      data%hel_state (:,2) = [1,-1, 1,-1, 0]

      allocate (data%col_state (2, 5, 2))
      data%col_state (:,:,1) = &
           reshape ([[1, 0], [2,-1], [3, 0], [2,-3], [0,0]], [2,5])
      data%col_state (:,:,2) = &
           reshape ([[1, 0], [2,-3], [3, 0], [2,-1], [0,0]], [2,5])

      allocate (data%ghost_flag (5, 2))
      data%ghost_flag(1:4,:) = .false.
      data%ghost_flag(5,:) = .true.

    end associate

    write (u, "(A)")  "* Set up the interaction"
    write (u, "(A)")

    call reset_interaction_counter ()
    call term%setup_interaction (core, model)
    call term%int%basic_write (u)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_13"
  end subroutine processes_13

  subroutine processes_14 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable, target :: process
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    type(sf_config_t), dimension(:), allocatable :: sf_config
    type(sf_channel_t), dimension(3) :: sf_channel

    write (u, "(A)")  "* Test output: processes_14"
    write (u, "(A)")  "*   Purpose: initialize a process with &
         &structure functions"
    write (u, "(A)")  "*            and compute MD5 sum"
    write (u, "(A)")

    write (u, "(A)")  "* Build and initialize a process object"
    write (u, "(A)")

    libname = "processes7"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib)
    call lib%compute_md5sum ()

    call model%init_test ()

    allocate (process)
    call process%init (procname, lib, os_data, model)

    call process%setup_test_cores ()
    allocate (phs_test_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    write (u, "(A)")  "* Set beam, structure functions, and mappings"
    write (u, "(A)")

    sqrts = 1000
    call process%setup_beams_sqrts (sqrts, i_core = 1)
    call process%configure_phs ()

    pdg_in = 25
    allocate (sf_test_data_t :: data)
    select type (data)
    type is (sf_test_data_t)
       call data%init (process%get_model_ptr (), pdg_in)
    end select

    call process%test_allocate_sf_channels (3)

    allocate (sf_config (2))
    call sf_config(1)%init ([1], data)
    call sf_config(2)%init ([2], data)
    call process%init_sf_chain (sf_config)
    deallocate (sf_config)

    call sf_channel(1)%init (2)
    call process%set_sf_channel (1, sf_channel(1))

    call sf_channel(2)%init (2)
    call sf_channel(2)%activate_mapping ([1,2])
    call process%set_sf_channel (2, sf_channel(2))

    call sf_channel(3)%init (2)
    call sf_channel(3)%set_s_mapping ([1,2])
    call process%set_sf_channel (3, sf_channel(3))

    call process%setup_mci (dispatch_mci_empty)

    call process%compute_md5sum ()

    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_14"

  end subroutine processes_14

  subroutine processes_15 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable, target :: process
    class(phs_config_t), allocatable :: phs_config_template
    type(process_instance_t), allocatable, target :: process_instance
    type(particle_set_t) :: pset

    write (u, "(A)")  "* Test output: processes_15"
    write (u, "(A)")  "*   Purpose: initialize a decay process object"
    write (u, "(A)")

    write (u, "(A)")  "* Build and load a test library with one process"
    write (u, "(A)")

    libname = "processes15"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib, scattering = .false., &
         decay = .true.)

    call model%init_test ()
    call model%set_par (var_str ("ff"), 0.4_default)
    call model%set_par (var_str ("mf"), &
         model%get_real (var_str ("ff")) * model%get_real (var_str ("ms")))

    write (u, "(A)")  "* Initialize a process object"
    write (u, "(A)")

    allocate (process)
    call process%init (procname, lib, os_data, model)

    call process%setup_test_cores ()
    allocate (phs_single_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    write (u, "(A)")  "* Prepare a trivial beam setup"
    write (u, "(A)")

    call process%setup_beams_decay (i_core = 1)
    call process%configure_phs ()
    call process%setup_mci (dispatch_mci_empty)

    write (u, "(A)")  "* Complete process initialization"
    write (u, "(A)")

    call process%setup_terms ()
    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Create a process instance"
    write (u, "(A)")

    call reset_interaction_counter (3)

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Inject a set of random numbers"
    write (u, "(A)")

    call process_instance%choose_mci (1)
    call process_instance%set_mcpar ([0._default, 0._default])
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Set up hard kinematics"
    write (u, "(A)")

    call process_instance%select_channel (1)
    call process_instance%compute_seed_kinematics ()
    call process_instance%compute_hard_kinematics ()

    write (u, "(A)")  "* Evaluate matrix element and square"
    write (u, "(A)")

    call process_instance%compute_eff_kinematics ()
    call process_instance%evaluate_expressions ()
    call process_instance%compute_other_channels ()
    call process_instance%evaluate_trace ()
    call process_instance%write (u)

    call process_instance%get_trace (pset, 1)
    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Particle content:"
    write (u, "(A)")

    call write_separator (u)
    call pset%write (u)
    call write_separator (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover process instance"
    write (u, "(A)")

    call reset_interaction_counter (3)

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%choose_mci (1)
    call process_instance%set_trace (pset, 1, check_match = .false.)
    call process_instance%recover (1, 1, .true., .true.)
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call pset%final ()
    call process_instance%final ()
    deallocate (process_instance)

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_15"

  end subroutine processes_15

  subroutine processes_16 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable, target :: process
    class(phs_config_t), allocatable :: phs_config_template
    type(process_instance_t), allocatable, target :: process_instance

    write (u, "(A)")  "* Test output: processes_16"
    write (u, "(A)")  "*   Purpose: integrate a process without &
         &structure functions"
    write (u, "(A)")  "*            in a multi-channel configuration"
    write (u, "(A)")

    write (u, "(A)")  "* Build and initialize a process object"
    write (u, "(A)")

    libname = "processes16"
    procname = libname

    call os_data%init ()
    call prc_test_create_library (libname, lib, scattering = .false., &
         decay = .true.)

    call reset_interaction_counter ()

    call model%init_test ()
    call model%set_par (var_str ("ff"), 0.4_default)
    call model%set_par (var_str ("mf"), &
         model%get_real (var_str ("ff")) * model%get_real (var_str ("ms")))

    allocate (process)
    call process%init (procname, lib, os_data, model)

    call process%setup_test_cores ()
    allocate (phs_single_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    write (u, "(A)")  "* Prepare a trivial beam setup"
    write (u, "(A)")

    call process%setup_beams_decay (i_core = 1)
    call process%configure_phs ()

    call process%setup_mci (dispatch_mci_test_midpoint)

    write (u, "(A)")  "* Complete process initialization"
    write (u, "(A)")

    call process%setup_terms ()
    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Create a process instance"
    write (u, "(A)")

    allocate (process_instance)
    call process_instance%init (process)

    write (u, "(A)")  "* Integrate with default test parameters"
    write (u, "(A)")

    call process_instance%integrate (1, n_it=1, n_calls=10000)
    call process%final_integration (1)

    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A,ES13.7)")  " Integral divided by phs factor = ", &
         process%get_integral (1) &
         / process_instance%kin(1)%phs_factor

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call process_instance%final ()
    deallocate (process_instance)

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_16"

  end subroutine processes_16

  subroutine processes_17 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_t), allocatable, target :: process
    class(phs_config_t), allocatable :: phs_config_template
    type(process_instance_t), allocatable, target :: process_instance
    type(particle_set_t) :: pset
    type(flavor_t) :: flv_beam
    real(default) :: m, p, E

    write (u, "(A)")  "* Test output: processes_17"
    write (u, "(A)")  "*   Purpose: initialize a decay process object"
    write (u, "(A)")

    write (u, "(A)")  "* Build and load a test library with one process"
    write (u, "(A)")

    libname = "processes17"
    procname = libname

    call os_data%init ()

    call prc_test_create_library (libname, lib, scattering = .false., &
         decay = .true.)

    write (u, "(A)")  "* Initialize a process object"
    write (u, "(A)")

    call model%init_test ()
    call model%set_par (var_str ("ff"), 0.4_default)
    call model%set_par (var_str ("mf"), &
         model%get_real (var_str ("ff")) * model%get_real (var_str ("ms")))

    allocate (process)
    call process%init (procname, lib, os_data, model)

    call process%setup_test_cores ()
    allocate (phs_single_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    write (u, "(A)")  "* Prepare a trivial beam setup"
    write (u, "(A)")

    call process%setup_beams_decay (rest_frame = .false., i_core = 1)
    call process%configure_phs ()
    call process%setup_mci (dispatch_mci_empty)

    write (u, "(A)")  "* Complete process initialization"
    write (u, "(A)")

    call process%setup_terms ()
    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Create a process instance"
    write (u, "(A)")

    call reset_interaction_counter (3)

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Set parent momentum and random numbers"
    write (u, "(A)")

    call process_instance%choose_mci (1)
    call process_instance%set_mcpar ([0._default, 0._default])

    call flv_beam%init (25, process%get_model_ptr ())
    m = flv_beam%get_mass ()
    p = 3 * m / 4
    E = sqrt (m**2 + p**2)
    call process_instance%set_beam_momenta ([vector4_moving (E, p, 3)])

    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Set up hard kinematics"
    write (u, "(A)")

    call process_instance%select_channel (1)
    call process_instance%compute_seed_kinematics ()
    call process_instance%compute_hard_kinematics ()

    write (u, "(A)")  "* Evaluate matrix element and square"
    write (u, "(A)")

    call process_instance%compute_eff_kinematics ()
    call process_instance%evaluate_expressions ()
    call process_instance%compute_other_channels ()
    call process_instance%evaluate_trace ()
    call process_instance%write (u)

    call process_instance%get_trace (pset, 1)
    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Particle content:"
    write (u, "(A)")

    call write_separator (u)
    call pset%write (u)
    call write_separator (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recover process instance"
    write (u, "(A)")

    call reset_interaction_counter (3)

    allocate (process_instance)
    call process_instance%init (process)

    call process_instance%choose_mci (1)
    call process_instance%set_trace (pset, 1, check_match = .false.)
    call process_instance%recover (1, 1, .true., .true.)
    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call pset%final ()
    call process_instance%final ()
    deallocate (process_instance)

    call process%final ()
    deallocate (process)

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_17"

  end subroutine processes_17

  subroutine processes_18 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(string_t) :: model_name
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    class(vars_t), pointer :: vars
    type(process_t), pointer :: process
    type(resonance_history_set_t) :: res_set
    integer :: i

    write (u, "(A)")  "* Test output: processes_18"
    write (u, "(A)")  "*   Purpose: extra resonance histories"
    write (u, "(A)")

    write (u, "(A)")  "* Build and load a test library with one process"
    write (u, "(A)")

    libname = "processes_18_lib"
    procname = "processes_18_p"

    call os_data%init ()

    call syntax_phs_forest_init ()

    model_name = "SM"
    model => null ()
    call prepare_model (model, model_name, vars)

    write (u, "(A)")  "* Initialize a process library with one process"
    write (u, "(A)")

    select type (model)
    class is (model_t)
       call prepare_resonance_test_library (lib, libname, procname, model, os_data, u)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Initialize a process object with phase space"

    allocate (process)
    select type (model)
    class is (model_t)
       call prepare_resonance_test_process (process, lib, procname, model, os_data)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Extract resonance history set"
    write (u, "(A)")

    call process%extract_resonance_history_set (res_set)
    call res_set%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call process%final ()
    deallocate (process)

    call model%final ()
    deallocate (model)

    call syntax_phs_forest_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_18"

  end subroutine processes_18

  subroutine processes_19 (u)
    integer, intent(in) :: u
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    type(process_t), allocatable, target :: process
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts
    type(process_instance_t) :: process_instance
    class(process_instance_hook_t), allocatable, target :: process_instance_hook, process_instance_hook2
    type(particle_set_t) :: pset

    write (u, "(A)")  "* Test output: processes_19"
    write (u, "(A)")  "*   Purpose: allocate process instance &
         &and add an after evaluate hook"
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)")  "* Allocate a process instance"
    write (u, "(A)")

    call process_instance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate hook and add to process instance"
    write (u, "(A)")

    allocate (process_instance_hook_test_t :: process_instance_hook)
    call process_instance%append_after_hook (process_instance_hook)

    allocate (process_instance_hook_test_t :: process_instance_hook2)
    call process_instance%append_after_hook (process_instance_hook2)

    select type (process_instance_hook)
    type is (process_instance_hook_test_t)
       process_instance_hook%unit = u
       process_instance_hook%name = "Hook 1"
    end select
    select type (process_instance_hook2)
    type is (process_instance_hook_test_t)
       process_instance_hook2%unit = u
       process_instance_hook2%name = "Hook 2"
    end select

    write (u, "(A)")  "* Evaluate matrix element and square"
    write (u, "(A)")

    call process_instance%evaluate_after_hook ()

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call process_instance_hook%final ()
    deallocate (process_instance_hook)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: processes_19"

  end subroutine processes_19


  subroutine dispatch_mci_empty (mci, var_list, process_id, is_nlo)
    class(mci_t), allocatable, intent(out) :: mci
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_id
    logical, intent(in), optional :: is_nlo
  end subroutine dispatch_mci_empty

  subroutine dispatch_mci_test3 (mci, var_list, process_id, is_nlo)
    class(mci_t), allocatable, intent(out) :: mci
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_id
    logical, intent(in), optional :: is_nlo
    allocate (mci_test_t :: mci)
    select type (mci)
    type is (mci_test_t)
       call mci%set_dimensions (2, 2)
       call mci%set_divisions (100)
    end select
  end subroutine dispatch_mci_test3

  subroutine dispatch_mci_test10 (mci, var_list, process_id, is_nlo)
    class(mci_t), allocatable, intent(out) :: mci
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_id
    logical, intent(in), optional :: is_nlo
    allocate (mci_test_t :: mci)
    select type (mci)
    type is (mci_test_t);  call mci%set_divisions (100)
    end select
  end subroutine dispatch_mci_test10

  subroutine prepare_test_process &
       (process, process_instance, model, var_list, run_id)
    type(process_t), intent(out), target :: process
    type(process_instance_t), intent(out), target :: process_instance
    class(model_data_t), intent(in), target :: model
    type(var_list_t), intent(inout), optional :: var_list
    type(string_t), intent(in), optional :: run_id
    type(process_library_t), target :: lib
    type(string_t) :: libname
    type(string_t) :: procname
    type(os_data_t) :: os_data
    type(model_t), allocatable, target :: process_model
    class(mci_t), pointer :: mci
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts
    libname = "processes_test"
    procname = libname
    call os_data%init ()
    call prc_test_create_library (libname, lib)
    call reset_interaction_counter ()
    allocate (process_model)
    call process_model%init (model%get_name (), &
         model%get_n_real (), &
         model%get_n_complex (), &
         model%get_n_field (), &
         model%get_n_vtx ())
    call process_model%copy_from (model)
    call process%init (procname, lib, os_data, process_model, var_list)
    if (present (run_id))  call process%set_run_id (run_id)
    call process%setup_test_cores ()
    allocate (phs_test_config_t :: phs_config_template)
    call process%init_components (phs_config_template)
    sqrts = 1000
    call process%setup_beams_sqrts (sqrts, i_core = 1)
    call process%configure_phs ()
    call process%setup_mci (dispatch_mci_test10)
    call process%setup_terms ()
    call process_instance%init (process)
    call process%test_get_mci_ptr (mci)
    select type (mci)
    type is (mci_test_t)
       ! This ensures that the next 'random' numbers are 0.3, 0.5, 0.7
       call mci%rng%init (3)
       ! Include the constant PHS factor in the stored maximum of the integrand
       call mci%set_max_factor (conv * twopi4 &
            / (2 * sqrt (lambda (sqrts **2, 125._default**2, 125._default**2))))
    end select
    call process%reset_library_ptr ()  ! avoid dangling pointer
    call process_model%final ()
  end subroutine prepare_test_process

  subroutine cleanup_test_process (process, process_instance)
    type(process_t), intent(inout) :: process
    type(process_instance_t), intent(inout) :: process_instance
    call process_instance%final ()
    call process%final ()
  end subroutine cleanup_test_process

  subroutine dispatch_mci_test_midpoint (mci, var_list, process_id, is_nlo)
    class(mci_t), allocatable, intent(out) :: mci
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_id
    logical, intent(in), optional :: is_nlo
    allocate (mci_midpoint_t :: mci)
  end subroutine dispatch_mci_test_midpoint

  subroutine prepare_resonance_test_library &
       (lib, libname, procname, model, os_data, u)
    type(process_library_t), target, intent(out) :: lib
    type(string_t), intent(in) :: libname
    type(string_t), intent(in) :: procname
    type(model_t), intent(in), target :: model
    type(os_data_t), intent(in) :: os_data
    integer, intent(in) :: u
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry

    call lib%init (libname)

    allocate (prt_in (2), prt_out (3))
    prt_in = [var_str ("e+"), var_str ("e-")]
    prt_out = [var_str ("d"), var_str ("ubar"), var_str ("W+")]

    allocate (template_me_def_t :: def)
    select type (def)
    type is (template_me_def_t)
       call def%init (model, prt_in, prt_out, unity = .false.)
    end select
    allocate (entry)
    call entry%init (procname, &
         model_name = model%get_name (), &
         n_in = 2, n_components = 1)
    call entry%import_component (1, n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("template"), &
         variant = def)
    call entry%write (u)

    call lib%append (entry)

    call lib%configure (os_data)
    call lib%write_makefile (os_data, force = .true., verbose = .false.)
    call lib%clean (os_data, distclean = .false.)
    call lib%write_driver (force = .true.)
    call lib%load (os_data)

  end subroutine prepare_resonance_test_library

  subroutine prepare_resonance_test_process &
       (process, lib, procname, model, os_data)
    class(process_t), intent(out), target :: process
    type(process_library_t), intent(in), target :: lib
    type(string_t), intent(in) :: procname
    type(model_t), intent(in), target :: model
    type(os_data_t), intent(in) :: os_data
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts

    call process%init (procname, lib, os_data, model)

    allocate (phs_wood_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    call process%setup_test_cores (type_string = var_str ("template"))

    sqrts = 1000
    call process%setup_beams_sqrts (sqrts, i_core = 1)
    call process%configure_phs ()
    call process%setup_mci (dispatch_mci_none)

    call process%setup_terms ()

  end subroutine prepare_resonance_test_process

  subroutine dispatch_mci_none (mci, var_list, process_id, is_nlo)
    class(mci_t), allocatable, intent(out) :: mci
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_id
    logical, intent(in), optional :: is_nlo
    allocate (mci_none_t :: mci)
  end subroutine dispatch_mci_none

  subroutine process_instance_hook_test_init (hook, var_list, instance, pdf_data)
    class(process_instance_hook_test_t), intent(inout), target :: hook
    type(var_list_t), intent(in) :: var_list
    class(process_instance_t), intent(in), target :: instance
    type(pdf_data_t), intent(in), optional :: pdf_data
  end subroutine process_instance_hook_test_init

  subroutine process_instance_hook_test_final (hook)
    class(process_instance_hook_test_t), intent(inout) :: hook
  end subroutine process_instance_hook_test_final

  subroutine process_instance_hook_test_evaluate (hook, instance)
    class(process_instance_hook_test_t), intent(inout) :: hook
    class(process_instance_t), intent(in), target :: instance
    write (hook%unit, "(A)") "Execute hook:"
    write (hook%unit, "(2X,A,1X,A,I0,A)") hook%name, "(", len (trim (hook%name)), ")"
  end subroutine process_instance_hook_test_evaluate


end module processes_uti

