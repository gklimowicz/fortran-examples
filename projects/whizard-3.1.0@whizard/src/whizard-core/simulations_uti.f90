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

module simulations_uti

    use kinds, only: default
    use kinds, only: i64
    use iso_varying_string, string_t => varying_string
    use io_units
    use format_defs, only: FMT_10, FMT_12
    use ifiles
    use lexers
    use parser
    use lorentz
    use flavors
    use interactions, only: reset_interaction_counter
    use process_libraries, only: process_library_t
    use prclib_stacks
    use phs_forests
    use event_base, only: generic_event_t
    use event_base, only: event_callback_t
    use particles, only: particle_set_t
    use eio_data
    use eio_base
    use eio_direct, only: eio_direct_t
    use eio_raw
    use eio_ascii
    use eio_dump
    use eio_callback
    use eval_trees
    use model_data, only: model_data_t
    use models
    use rt_data
    use event_streams
    use decays_ut, only: prepare_testbed
    use process, only: process_t
    use process_stacks, only: process_entry_t
    use process_configurations_ut, only: prepare_test_library
    use compilations, only: compile_library
    use integrations, only: integrate_process

    use simulations

    use restricted_subprocesses_uti, only: prepare_resonance_test_library

  implicit none
  private

  public :: simulations_1
  public :: simulations_2
  public :: simulations_3
  public :: simulations_4
  public :: simulations_5
  public :: simulations_6
  public :: simulations_7
  public :: simulations_8
  public :: simulations_9
  public :: simulations_10
  public :: simulations_11
  public :: simulations_12
  public :: simulations_13
  public :: simulations_14
  public :: simulations_15

  type, extends (event_callback_t) :: simulations_13_callback_t
     integer :: u
   contains
     procedure :: write => simulations_13_callback_write
     procedure :: proc => simulations_13_callback
  end type simulations_13_callback_t


contains

  subroutine simulations_1 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1, procname2
    type(rt_data_t), target :: global
    type(simulation_t), target :: simulation

    write (u, "(A)")  "* Test output: simulations_1"
    write (u, "(A)")  "*   Purpose: initialize simulation"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize processes"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_1a"
    procname1 = "simulation_1p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("simulations1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    procname2 = "sim_extra"

    call prepare_test_library (global, libname, 1, [procname2])
    call compile_library (libname, global)
    call global%set_string (var_str ("$run_id"), &
         var_str ("simulations2"), is_known = .true.)


    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call global%set_string (var_str ("$sample"), &
         var_str ("sim1"), is_known = .true.)
    call integrate_process (procname2, global, local_stack=.true.)

    call simulation%init ([procname1, procname2], .false., .true., global)
    call simulation%init_process_selector ()
    call simulation%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write the event record for the first process"
    write (u, "(A)")

    call simulation%write_event (u, i_prc = 1)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_1"

  end subroutine simulations_1

  subroutine simulations_2 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1
    type(rt_data_t), target :: global
    type(simulation_t), target :: simulation
    type(event_sample_data_t) :: data

    write (u, "(A)")  "* Test output: simulations_2"
    write (u, "(A)")  "*   Purpose: generate events for a single process"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize processes"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_2a"
    procname1 = "simulation_2p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("simulations1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call global%set_log (var_str ("?unweighted"), &
         .false., is_known = .true.)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    data = simulation%get_data ()
    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate three events"
    write (u, "(A)")

    call simulation%set_n_events_requested (3)
    call simulation%generate ()
    call simulation%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write the event record for the last event"
    write (u, "(A)")

    call simulation%write_event (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_2"

  end subroutine simulations_2

  subroutine simulations_3 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1
    type(rt_data_t), target :: global
    type(simulation_t), target :: simulation
    type(event_sample_data_t) :: data

    write (u, "(A)")  "* Test output: simulations_3"
    write (u, "(A)")  "*   Purpose: generate unweighted events &
         &for a single process"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize processes"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_3a"
    procname1 = "simulation_3p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("simulations1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    data = simulation%get_data ()
    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate three events"
    write (u, "(A)")

    call simulation%set_n_events_requested (3)
    call simulation%generate ()
    call simulation%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write the event record for the last event"
    write (u, "(A)")

    call simulation%write_event (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_3"

  end subroutine simulations_3

  subroutine simulations_4 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1
    type(rt_data_t), target :: global
    type(flavor_t) :: flv
    type(string_t) :: name
    type(simulation_t), target :: simulation
    type(event_sample_data_t) :: data

    write (u, "(A)")  "* Test output: simulations_4"
    write (u, "(A)")  "*   Purpose: generate events for a single process &
         &with structure functions"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize processes"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_4a"
    procname1 = "simulation_4p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("wood"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .true., is_known = .true.)
    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%model_set_real (var_str ("ms"), &
         0._default)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call reset_interaction_counter ()

    call flv%init (25, global%model)
    name = flv%get_name ()

    call global%beam_structure%init_sf ([name, name], [1])
    call global%beam_structure%set_sf (1, 1, var_str ("sf_test_1"))

    write (u, "(A)")  "* Integrate"
    write (u, "(A)")

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call global%set_log (var_str ("?unweighted"), &
         .false., is_known = .true.)
    call global%set_string (var_str ("$sample"), &
         var_str ("simulations4"), is_known = .true.)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    data = simulation%get_data ()
    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate three events"
    write (u, "(A)")

    call simulation%set_n_events_requested (3)
    call simulation%generate ()
    call simulation%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write the event record for the last event"
    write (u, "(A)")

    call simulation%write_event (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_4"

  end subroutine simulations_4

  subroutine simulations_5 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1, sample
    type(rt_data_t), target :: global
    class(eio_t), allocatable :: eio
    type(simulation_t), allocatable, target :: simulation

    write (u, "(A)")  "* Test output: simulations_5"
    write (u, "(A)")  "*   Purpose: generate events for a single process"
    write (u, "(A)")  "*            write to file and reread"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize processes"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_5a"
    procname1 = "simulation_5p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("simulations5"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call global%set_log (var_str ("?unweighted"), &
         .false., is_known = .true.)
    sample = "simulations5"
    call global%set_string (var_str ("$sample"), &
         sample, is_known = .true.)
    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    write (u, "(A)")  "* Initialize raw event file"
    write (u, "(A)")

    allocate (eio_raw_t :: eio)
    call eio%init_out (sample)

    write (u, "(A)")  "* Generate an event"
    write (u, "(A)")

    call simulation%set_n_events_requested (1)
    call simulation%generate ()
    call simulation%write_event (u)
    call simulation%write_event (eio)

    call eio%final ()
    deallocate (eio)
    call simulation%final ()
    deallocate (simulation)

    write (u, "(A)")
    write (u, "(A)")  "* Re-read the event from file"
    write (u, "(A)")

    call global%set_log (var_str ("?update_sqme"), &
         .true., is_known = .true.)
    call global%set_log (var_str ("?update_weight"), &
         .true., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()
    allocate (eio_raw_t :: eio)
    call eio%init_in (sample)

    call simulation%read_event (eio)
    call simulation%write_event (u)

    write (u, "(A)")
    write (u, "(A)")  "* Recalculate process instance"
    write (u, "(A)")

    call simulation%recalculate ()
    call simulation%evaluate_expressions ()
    call simulation%write_event (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio%final ()
    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_5"

  end subroutine simulations_5

  subroutine simulations_6 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1, sample
    type(rt_data_t), target :: global
    class(eio_t), allocatable :: eio
    type(simulation_t), allocatable, target :: simulation
    type(flavor_t) :: flv
    type(string_t) :: name

    write (u, "(A)")  "* Test output: simulations_6"
    write (u, "(A)")  "*   Purpose: generate events for a single process"
    write (u, "(A)")  "*            write to file and reread"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and integrate"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_6"
    procname1 = "simulation_6p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("wood"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .true., is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%model_set_real (var_str ("ms"), &
         0._default)

    call flv%init (25, global%model)
    name = flv%get_name ()

    call global%beam_structure%init_sf ([name, name], [1])
    call global%beam_structure%set_sf (1, 1, var_str ("sf_test_1"))

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call reset_interaction_counter ()

    call global%set_log (var_str ("?unweighted"), &
         .false., is_known = .true.)
    sample = "simulations6"
    call global%set_string (var_str ("$sample"), &
         sample, is_known = .true.)
    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    write (u, "(A)")  "* Initialize raw event file"
    write (u, "(A)")

    allocate (eio_raw_t :: eio)
    call eio%init_out (sample)

    write (u, "(A)")  "* Generate an event"
    write (u, "(A)")

    call simulation%set_n_events_requested (1)
    call simulation%generate ()
    call pacify (simulation)
    call simulation%write_event (u, verbose = .true., testflag = .true.)
    call simulation%write_event (eio)

    call eio%final ()
    deallocate (eio)
    call simulation%final ()
    deallocate (simulation)

    write (u, "(A)")
    write (u, "(A)")  "* Re-read the event from file"
    write (u, "(A)")

    call reset_interaction_counter ()

    call global%set_log (var_str ("?update_sqme"), &
         .true., is_known = .true.)
    call global%set_log (var_str ("?update_weight"), &
         .true., is_known = .true.)

    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()
    allocate (eio_raw_t :: eio)
    call eio%init_in (sample)

    call simulation%read_event (eio)
    call simulation%write_event (u, verbose = .true., testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Recalculate process instance"
    write (u, "(A)")

    call simulation%recalculate ()
    call simulation%evaluate_expressions ()
    call simulation%write_event (u, verbose = .true., testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio%final ()
    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_6"

  end subroutine simulations_6

  subroutine simulations_7 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1, sample
    type(rt_data_t), target :: global
    type(string_t), dimension(0) :: empty_string_array
    type(event_sample_data_t) :: data
    type(event_stream_array_t) :: es_array
    type(simulation_t), allocatable, target :: simulation
    type(flavor_t) :: flv
    type(string_t) :: name

    write (u, "(A)")  "* Test output: simulations_7"
    write (u, "(A)")  "*   Purpose: generate events for a single process"
    write (u, "(A)")  "*            write to file and reread"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and integrate"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_7"
    procname1 = "simulation_7p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("wood"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .true., is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%model_set_real (var_str ("ms"), &
         0._default)

    call flv%init (25, global%model)
    name = flv%get_name ()

    call global%beam_structure%init_sf ([name, name], [1])
    call global%beam_structure%set_sf (1, 1, var_str ("sf_test_1"))

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call reset_interaction_counter ()

    call global%set_log (var_str ("?unweighted"), &
         .false., is_known = .true.)
    sample = "simulations7"
    call global%set_string (var_str ("$sample"), &
         sample, is_known = .true.)
    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    write (u, "(A)")  "* Initialize raw event file"
    write (u, "(A)")

    data%md5sum_prc = simulation%get_md5sum_prc ()
    data%md5sum_cfg = simulation%get_md5sum_cfg ()
    call es_array%init (sample, [var_str ("raw")], global, data)

    write (u, "(A)")  "* Generate an event"
    write (u, "(A)")

    call simulation%set_n_events_requested (1)
    call simulation%generate (es_array)

    call es_array%final ()
    call simulation%final ()
    deallocate (simulation)

    write (u, "(A)")  "* Re-read the event from file and generate another one"
    write (u, "(A)")

    call global%set_log (&
         var_str ("?rebuild_events"), .false., is_known = .true.)

    call reset_interaction_counter ()

    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    data%md5sum_prc = simulation%get_md5sum_prc ()
    data%md5sum_cfg = simulation%get_md5sum_cfg ()
    call es_array%init (sample, empty_string_array, global, data, &
         input = var_str ("raw"))

    call simulation%set_n_events_requested (2)
    call simulation%generate (es_array)

    call pacify (simulation)
    call simulation%write_event (u, verbose = .true.)

    call es_array%final ()
    call simulation%final ()
    deallocate (simulation)


    write (u, "(A)")
    write (u, "(A)")  "* Re-read both events from file"
    write (u, "(A)")

    call reset_interaction_counter ()

    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    data%md5sum_prc = simulation%get_md5sum_prc ()
    data%md5sum_cfg = simulation%get_md5sum_cfg ()
    call es_array%init (sample, empty_string_array, global, data, &
         input = var_str ("raw"))

    call simulation%set_n_events_requested (2)
    call simulation%generate (es_array)

    call pacify (simulation)
    call simulation%write_event (u, verbose = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call es_array%final ()
    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_7"

  end subroutine simulations_7

  subroutine simulations_8 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1, sample
    type(rt_data_t), target :: global
    type(string_t), dimension(0) :: empty_string_array
    type(event_sample_data_t) :: data
    type(event_stream_array_t) :: es_array
    type(simulation_t), allocatable, target :: simulation
    type(flavor_t) :: flv
    type(string_t) :: name

    write (u, "(A)")  "* Test output: simulations_8"
    write (u, "(A)")  "*   Purpose: generate events for a single process"
    write (u, "(A)")  "*            write to file and rescan"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and integrate"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_8"
    procname1 = "simulation_8p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("wood"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .true., is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%model_set_real (var_str ("ms"), &
         0._default)

    call flv%init (25, global%model)
    name = flv%get_name ()

    call global%beam_structure%init_sf ([name, name], [1])
    call global%beam_structure%set_sf (1, 1, var_str ("sf_test_1"))

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call reset_interaction_counter ()

    call global%set_log (var_str ("?unweighted"), &
         .false., is_known = .true.)
    sample = "simulations8"
    call global%set_string (var_str ("$sample"), &
         sample, is_known = .true.)
    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    write (u, "(A)")  "* Initialize raw event file"
    write (u, "(A)")

    data%md5sum_prc = simulation%get_md5sum_prc ()
    data%md5sum_cfg = simulation%get_md5sum_cfg ()
    write (u, "(1x,A,A,A)")  "MD5 sum (proc)   = '", data%md5sum_prc, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (config) = '", data%md5sum_cfg, "'"
    call es_array%init (sample, [var_str ("raw")], global, &
         data)

    write (u, "(A)")
    write (u, "(A)")  "* Generate an event"
    write (u, "(A)")

    call simulation%set_n_events_requested (1)
    call simulation%generate (es_array)

    call pacify (simulation)
    call simulation%write_event (u, verbose = .true., testflag = .true.)

    call es_array%final ()
    call simulation%final ()
    deallocate (simulation)

    write (u, "(A)")
    write (u, "(A)")  "* Re-read the event from file"
    write (u, "(A)")

    call reset_interaction_counter ()

    allocate (simulation)
    call simulation%init ([procname1], .false., .false., global)
    call simulation%init_process_selector ()

    data%md5sum_prc = simulation%get_md5sum_prc ()
    data%md5sum_cfg = ""
    write (u, "(1x,A,A,A)")  "MD5 sum (proc)   = '", data%md5sum_prc, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (config) = '", data%md5sum_cfg, "'"
    call es_array%init (sample, empty_string_array, global, data, &
         input = var_str ("raw"), input_sample = sample, allow_switch = .false.)

    call simulation%rescan (1, es_array, global = global)

    write (u, "(A)")

    call pacify (simulation)
    call simulation%write_event (u, verbose = .true., testflag = .true.)

    call es_array%final ()
    call simulation%final ()
    deallocate (simulation)

    write (u, "(A)")
    write (u, "(A)")  "* Re-read again and recalculate"
    write (u, "(A)")

    call reset_interaction_counter ()

    call global%set_log (var_str ("?update_sqme"), &
         .true., is_known = .true.)
    call global%set_log (var_str ("?update_event"), &
         .true., is_known = .true.)

    allocate (simulation)
    call simulation%init ([procname1], .false., .false., global)
    call simulation%init_process_selector ()

    data%md5sum_prc = simulation%get_md5sum_prc ()
    data%md5sum_cfg = ""
    write (u, "(1x,A,A,A)")  "MD5 sum (proc)   = '", data%md5sum_prc, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (config) = '", data%md5sum_cfg, "'"
    call es_array%init (sample, empty_string_array, global, data, &
         input = var_str ("raw"), input_sample = sample, allow_switch = .false.)

    call simulation%rescan (1, es_array, global = global)

    write (u, "(A)")

    call pacify (simulation)
    call simulation%write_event (u, verbose = .true., testflag = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call es_array%final ()
    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_8"

  end subroutine simulations_8

  subroutine simulations_9 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1, sample
    type(rt_data_t), target :: global
    type(string_t), dimension(0) :: empty_string_array
    type(event_sample_data_t) :: data
    type(event_stream_array_t) :: es_array
    type(simulation_t), allocatable, target :: simulation
    type(flavor_t) :: flv
    type(string_t) :: name
    logical :: error

    write (u, "(A)")  "* Test output: simulations_9"
    write (u, "(A)")  "*   Purpose: generate events for a single process"
    write (u, "(A)")  "*            write to file and rescan"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and integrate"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))

    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_9"
    procname1 = "simulation_9p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("wood"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("vamp"), is_known = .true.)
    call global%set_log (var_str ("?use_vamp_equivalences"),&
         .true., is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%model_set_real (var_str ("ms"), &
         0._default)

    call flv%init (25, global%model)
    name = flv%get_name ()

    call global%beam_structure%init_sf ([name, name], [1])
    call global%beam_structure%set_sf (1, 1, var_str ("sf_test_1"))

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call reset_interaction_counter ()

    call global%set_log (var_str ("?unweighted"), &
         .false., is_known = .true.)
    sample = "simulations9"
    call global%set_string (var_str ("$sample"), &
         sample, is_known = .true.)
    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    call simulation%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize raw event file"
    write (u, "(A)")

    data%md5sum_prc = simulation%get_md5sum_prc ()
    data%md5sum_cfg = simulation%get_md5sum_cfg ()
    write (u, "(1x,A,A,A)")  "MD5 sum (proc)   = '", data%md5sum_prc, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (config) = '", data%md5sum_cfg, "'"
    call es_array%init (sample, [var_str ("raw")], global, &
         data)

    write (u, "(A)")
    write (u, "(A)")  "* Generate an event"
    write (u, "(A)")

    call simulation%set_n_events_requested (1)
    call simulation%generate (es_array)

    call es_array%final ()
    call simulation%final ()
    deallocate (simulation)

    write (u, "(A)")  "* Initialize event generation for different parameters"
    write (u, "(A)")

    call reset_interaction_counter ()

    allocate (simulation)
    call simulation%init ([procname1, procname1], .false., .false., global)
    call simulation%init_process_selector ()

    call simulation%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Attempt to re-read the events (should fail)"
    write (u, "(A)")

    data%md5sum_prc = simulation%get_md5sum_prc ()
    data%md5sum_cfg = ""
    write (u, "(1x,A,A,A)")  "MD5 sum (proc)   = '", data%md5sum_prc, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (config) = '", data%md5sum_cfg, "'"
    call es_array%init (sample, empty_string_array, global, data, &
         input = var_str ("raw"), input_sample = sample, &
         allow_switch = .false., error = error)

    write (u, "(1x,A,L1)")  "error = ", error

    call simulation%rescan (1, es_array, global = global)

    call es_array%final ()
    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_9"

  end subroutine simulations_9

  subroutine simulations_10 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1, expr_text
    type(rt_data_t), target :: global
    type(rt_data_t), dimension(1), target :: alt_env
    type(ifile_t) :: ifile
    type(stream_t) :: stream
    type(parse_tree_t) :: pt_weight
    type(simulation_t), target :: simulation
    type(event_sample_data_t) :: data

    write (u, "(A)")  "* Test output: simulations_10"
    write (u, "(A)")  "*   Purpose: reweight event"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize processes"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_pexpr_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_10a"
    procname1 = "simulation_10p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("simulations1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize alternative environment with custom weight"
    write (u, "(A)")

    call alt_env(1)%local_init (global)
    call alt_env(1)%activate ()

    expr_text = "2"
    write (u, "(A,A)")  "weight = ", char (expr_text)
    write (u, *)

    call ifile_clear (ifile)
    call ifile_append (ifile, expr_text)
    call stream_init (stream, ifile)
    call parse_tree_init_expr (pt_weight, stream, .true.)
    call stream_final (stream)
    alt_env(1)%pn%weight_expr => pt_weight%get_root_ptr ()
    call alt_env(1)%write_expr (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call global%set_log (var_str ("?unweighted"), &
         .false., is_known = .true.)
    call simulation%init ([procname1], .true., .true., global, alt_env=alt_env)
    call simulation%init_process_selector ()

    data = simulation%get_data ()
    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate an event"
    write (u, "(A)")

    call simulation%set_n_events_requested (1)
    call simulation%generate ()
    call simulation%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write the event record for the last event"
    write (u, "(A)")

    call simulation%write_event (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write the event record for the alternative setup"
    write (u, "(A)")

    call simulation%write_alt_event (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call simulation%final ()
    call global%final ()

    call syntax_model_file_final ()
    call syntax_pexpr_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_10"

  end subroutine simulations_10

  subroutine simulations_11 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global
    type(prclib_entry_t), pointer :: lib
    type(string_t) :: prefix, procname1, procname2
    type(simulation_t), target :: simulation

    write (u, "(A)")  "* Test output: simulations_11"
    write (u, "(A)")  "*   Purpose: apply decay"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize processes"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    allocate (lib)
    call global%add_prclib (lib)

    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    prefix = "simulation_11"
    procname1 = prefix // "_p"
    procname2 = prefix // "_d"
    call prepare_testbed &
         (global%prclib, global%process_stack, &
         prefix, global%os_data, &
         scattering=.true., decay=.true.)

    call global%select_model (var_str ("Test"))
    call global%model%set_par (var_str ("ff"), 0.4_default)
    call global%model%set_par (var_str ("mf"), &
         global%model%get_real (var_str ("ff")) &
         * global%model%get_real (var_str ("ms")))
    call global%model%set_unstable (25, [procname2])

    write (u, "(A)")  "* Initialize simulation object"
    write (u, "(A)")

    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    write (u, "(A)")  "* Generate event"
    write (u, "(A)")

    call simulation%set_n_events_requested (1)
    call simulation%generate ()
    call simulation%write (u)

    write (u, *)

    call simulation%write_event (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"
    write (u, "(A)")

    call simulation%final ()
    call global%final ()

    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_11"

  end subroutine simulations_11

  subroutine simulations_12 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1, sample
    type(rt_data_t), target :: global
    class(eio_t), allocatable :: eio
    type(simulation_t), allocatable, target :: simulation
    type(flavor_t) :: flv
    integer :: i_evt

    write (u, "(A)")  "* Test output: simulations_12"
    write (u, "(A)")  "*   Purpose: generate events for a single process"
    write (u, "(A)")  "*            and write to split event files"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and integrate"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_12"
    procname1 = "simulation_12p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%model_set_real (var_str ("ms"), &
         0._default)

    call flv%init (25, global%model)

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call global%set_log (var_str ("?unweighted"), &
         .false., is_known = .true.)
    sample = "simulations_12"
    call global%set_string (var_str ("$sample"), &
         sample, is_known = .true.)
    call global%set_int (var_str ("sample_split_n_evt"), &
         2, is_known = .true.)
    call global%set_int (var_str ("sample_split_index"), &
         42, is_known = .true.)
    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    call simulation%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize ASCII event file"
    write (u, "(A)")

    allocate (eio_ascii_short_t :: eio)
    select type (eio)
    class is (eio_ascii_t);  call eio%set_parameters ()
    end select
    call eio%init_out (sample, data = simulation%get_data ())

    write (u, "(A)")  "* Generate 5 events, distributed among three files"

    do i_evt = 1, 5
       call simulation%set_n_events_requested (1)
       call simulation%generate ()
       call simulation%write_event (eio)
    end do

    call eio%final ()
    deallocate (eio)
    call simulation%final ()
    deallocate (simulation)

    write (u, *)
    call display_file ("simulations_12.42.short.evt", u)
    write (u, *)
    call display_file ("simulations_12.43.short.evt", u)
    write (u, *)
    call display_file ("simulations_12.44.short.evt", u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_12"

  end subroutine simulations_12

  subroutine simulations_13 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, procname1, sample
    type(rt_data_t), target :: global
    class(eio_t), allocatable :: eio
    type(simulation_t), allocatable, target :: simulation
    type(flavor_t) :: flv
    integer :: i_evt
    type(simulations_13_callback_t) :: event_callback

    write (u, "(A)")  "* Test output: simulations_13"
    write (u, "(A)")  "*   Purpose: generate events for a single process"
    write (u, "(A)")  "*            and execute callback"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize process and integrate"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)

    libname = "simulation_13"
    procname1 = "simulation_13p"

    call prepare_test_library (global, libname, 1, [procname1])
    call compile_library (libname, global)

    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call global%set_string (var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call global%set_string (var_str ("$integration_method"),&
         var_str ("midpoint"), is_known = .true.)
    call global%set_log (var_str ("?vis_history"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?integration_timer"),&
         .false., is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)

    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)

    call flv%init (25, global%model)

    call global%it_list%init ([1], [1000])

    call global%set_string (var_str ("$run_id"), &
         var_str ("r1"), is_known = .true.)
    call integrate_process (procname1, global, local_stack=.true.)

    write (u, "(A)")  "* Initialize event generation"
    write (u, "(A)")

    call global%set_log (var_str ("?unweighted"), &
         .false., is_known = .true.)
    sample = "simulations_13"
    call global%set_string (var_str ("$sample"), &
         sample, is_known = .true.)

    allocate (simulation)
    call simulation%init ([procname1], .true., .true., global)
    call simulation%init_process_selector ()

    write (u, "(A)")  "* Prepare callback object"
    write (u, "(A)")

    event_callback%u = u
    call global%set_event_callback (event_callback)

    write (u, "(A)")  "* Initialize callback I/O object"
    write (u, "(A)")

    allocate (eio_callback_t :: eio)
    select type (eio)
    class is (eio_callback_t)
       call eio%set_parameters (callback = event_callback, &
            count_interval = 3)
    end select
    call eio%init_out (sample, data = simulation%get_data ())

    write (u, "(A)")  "* Generate 7 events, with callback every 3 events"
    write (u, "(A)")

    do i_evt = 1, 7
       call simulation%set_n_events_requested (1)
       call simulation%generate ()
       call simulation%write_event (eio)
    end do

    call eio%final ()
    deallocate (eio)
    call simulation%final ()
    deallocate (simulation)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_13"

  end subroutine simulations_13

  subroutine simulations_14 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, libname_generated
    type(string_t) :: procname
    type(string_t) :: model_name
    type(rt_data_t), target :: global
    type(prclib_entry_t), pointer :: lib_entry
    type(process_library_t), pointer :: lib
    class(model_t), pointer :: model
    class(model_data_t), pointer :: model_data
    type(simulation_t), target :: simulation
    type(particle_set_t) :: pset
    type(eio_direct_t) :: eio_in
    type(eio_dump_t) :: eio_out
    real(default) :: sqrts, mw, pp
    real(default), dimension(3) :: p3
    type(vector4_t), dimension(:), allocatable :: p
    real(default), dimension(:), allocatable :: m
    integer :: u_verbose, i
    real(default) :: sqme_proc
    real(default), dimension(:), allocatable :: sqme
    real(default) :: on_shell_limit
    integer, dimension(:), allocatable :: i_array
    real(default), dimension(:), allocatable :: prob_array

    write (u, "(A)")  "* Test output: simulations_14"
    write (u, "(A)")  "*   Purpose: construct resonant subprocesses &
         &in the simulation object"
    write (u, "(A)")

    write (u, "(A)")  "* Build and load a test library with one process"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    libname = "simulations_14_lib"
    procname = "simulations_14_p"

    call global%global_init ()
    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)
    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)
    call global%set_log (var_str ("?update_sqme"), &
         .true., is_known = .true.)
    call global%set_log (var_str ("?update_weight"), &
         .true., is_known = .true.)
    call global%set_log (var_str ("?update_event"), &
         .true., is_known = .true.)

    model_name = "SM"
    call global%select_model (model_name)
    allocate (model)
    call model%init_instance (global%model)
    model_data => model

    write (u, "(A)")  "* Initialize process library and process"
    write (u, "(A)")

    allocate (lib_entry)
    call lib_entry%init (libname)
    lib => lib_entry%process_library_t
    call global%add_prclib (lib_entry)

    call prepare_resonance_test_library &
         (lib, libname, procname, model_data, global, u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize simulation object &
         &with resonant subprocesses"
    write (u, "(A)")

    call global%set_log (var_str ("?resonance_history"), &
         .true., is_known = .true.)
    call global%set_real (var_str ("resonance_on_shell_limit"), &
         10._default, is_known = .true.)

    call simulation%init ([procname], &
         integrate=.false., generate=.false., local=global)

    call simulation%write_resonant_subprocess_data (u, 1)

    write (u, "(A)")
    write (u, "(A)")  "* Resonant subprocesses: generated library"
    write (u, "(A)")

    libname_generated = procname // "_R"
    lib => global%prclib_stack%get_library_ptr (libname_generated)
    if (associated (lib))  call lib%write (u, libpath=.false.)

    write (u, "(A)")
    write (u, "(A)")  "* Generated process stack"
    write (u, "(A)")

    call global%process_stack%show (u)

    write (u, "(A)")
    write (u, "(A)")  "* Particle set"
    write (u, "(A)")

    pset = simulation%get_hard_particle_set (1)
    call pset%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize object for direct access"
    write (u, "(A)")

    call eio_in%init_direct &
         (n_beam = 0, n_in = 2, n_rem = 0, n_vir = 0, n_out = 3, &
         pdg = [-11, 11, 1, -2, 24], model=global%model)
    call eio_in%set_selection_indices (1, 1, 1, 1)

    sqrts = global%get_rval (var_str ("sqrts"))
    mw = 80._default   ! deliberately slightly different from true mw
    pp = sqrt (sqrts**2 - 4 * mw**2) / 2

    allocate (p (5), m (5))
    p(1) = vector4_moving (sqrts/2, sqrts/2, 3)
    m(1) = 0
    p(2) = vector4_moving (sqrts/2,-sqrts/2, 3)
    m(2) = 0
    p3(1) = pp/2
    p3(2) = mw/2
    p3(3) = 0
    p(3) = vector4_moving (sqrts/4, vector3_moving (p3))
    m(3) = 0
    p3(2) = -mw/2
    p(4) = vector4_moving (sqrts/4, vector3_moving (p3))
    m(4) = 0
    p(5) = vector4_moving (sqrts/2,-pp, 1)
    m(5) = mw
    call eio_in%set_momentum (p, m**2)

    call eio_in%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Transfer and show particle set"
    write (u, "(A)")

    call simulation%read_event (eio_in)
    pset = simulation%get_hard_particle_set (1)
    call pset%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* (Re)calculate matrix element"
    write (u, "(A)")

    call simulation%recalculate (recover_phs = .false.)
    call simulation%evaluate_transforms ()

    write (u, "(A)")  "* Show event with sqme"
    write (u, "(A)")

    call eio_out%set_parameters (unit = u, &
         weights = .true., pacify = .true., compressed = .true.)
    call eio_out%init_out (var_str (""))
    call simulation%write_event (eio_out)

    write (u, "(A)")
    write (u, "(A)")  "* Write event to separate file &
         &'simulations_14_event_verbose.log'"

    u_verbose = free_unit ()
    open (unit = u_verbose, file = "simulations_14_event_verbose.log", &
         status = "replace", action = "write")
    call simulation%write (u_verbose)
    write (u_verbose, *)
    call simulation%write_event (u_verbose, verbose =.true., testflag = .true.)
    close (u_verbose)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_14"

  end subroutine simulations_14

  subroutine simulations_15 (u)
    integer, intent(in) :: u
    type(string_t) :: libname, libname_generated
    type(string_t) :: procname
    type(string_t) :: model_name
    type(rt_data_t), target :: global
    type(prclib_entry_t), pointer :: lib_entry
    type(process_library_t), pointer :: lib
    class(model_t), pointer :: model
    class(model_data_t), pointer :: model_data
    type(simulation_t), target :: simulation
    real(default) :: sqrts
    type(eio_dump_t) :: eio_out
    integer :: u_verbose

    write (u, "(A)")  "* Test output: simulations_15"
    write (u, "(A)")  "*   Purpose: generate event with resonant subprocess"
    write (u, "(A)")

    write (u, "(A)")  "* Build and load a test library with one process"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    libname = "simulations_15_lib"
    procname = "simulations_15_p"

    call global%global_init ()
    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_grids"), .true., intrinsic = .true.)
    call global%append_log (&
         var_str ("?rebuild_events"), .true., intrinsic = .true.)
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)
    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%set_log (var_str ("?recover_beams"), &
         .false., is_known = .true.)
    call global%set_log (var_str ("?update_sqme"), &
         .true., is_known = .true.)
    call global%set_log (var_str ("?update_weight"), &
         .true., is_known = .true.)
    call global%set_log (var_str ("?update_event"), &
         .true., is_known = .true.)
    call global%set_log (var_str ("?resonance_history"), &
         .true., is_known = .true.)
    call global%set_real (var_str ("resonance_on_shell_limit"), &
         10._default, is_known = .true.)

    model_name = "SM"
    call global%select_model (model_name)
    allocate (model)
    call model%init_instance (global%model)
    model_data => model

    write (u, "(A)")  "* Initialize process library and process"
    write (u, "(A)")

    allocate (lib_entry)
    call lib_entry%init (libname)
    lib => lib_entry%process_library_t
    call global%add_prclib (lib_entry)

    call prepare_resonance_test_library &
         (lib, libname, procname, model_data, global, u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize simulation object &
         &with resonant subprocesses"
    write (u, "(A)")

    call global%it_list%init ([1], [1000])
    call simulation%init ([procname], &
         integrate=.true., generate=.true., local=global)

    call simulation%write_resonant_subprocess_data (u, 1)

    write (u, "(A)")
    write (u, "(A)")  "* Generate event"
    write (u, "(A)")

    call simulation%init_process_selector ()
    call simulation%set_n_events_requested (1)
    call simulation%generate ()

    call eio_out%set_parameters (unit = u, &
         weights = .true., pacify = .true., compressed = .true.)
    call eio_out%init_out (var_str (""))
    call simulation%write_event (eio_out)

    write (u, "(A)")
    write (u, "(A)")  "* Write event to separate file &
         &'simulations_15_event_verbose.log'"

    u_verbose = free_unit ()
    open (unit = u_verbose, file = "simulations_15_event_verbose.log", &
         status = "replace", action = "write")
    call simulation%write (u_verbose)
    write (u_verbose, *)
    call simulation%write_event (u_verbose, verbose =.true., testflag = .true.)
    close (u_verbose)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call simulation%final ()
    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: simulations_15"

  end subroutine simulations_15


  subroutine display_file (file, u)
    use io_units, only: free_unit
    character(*), intent(in) :: file
    integer, intent(in) :: u
    character(256) :: buffer
    integer :: u_file
    write (u, "(3A)")  "* Contents of file '", file, "':"
    write (u, *)
    u_file = free_unit ()
    open (u_file, file = file, action = "read", status = "old")
    do
       read (u_file, "(A)", end = 1)  buffer
       write (u, "(A)")  trim (buffer)
    end do
1   continue
  end subroutine display_file

  subroutine simulations_13_callback_write (event_callback, unit)
    class(simulations_13_callback_t), intent(in) :: event_callback
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Hello"
  end subroutine simulations_13_callback_write

  subroutine simulations_13_callback (event_callback, i, event)
    class(simulations_13_callback_t), intent(in) :: event_callback
    integer(i64), intent(in) :: i
    class(generic_event_t), intent(in) :: event
    write (event_callback%u, "(A,I0)")  "hello event #", i
  end subroutine simulations_13_callback


end module simulations_uti

