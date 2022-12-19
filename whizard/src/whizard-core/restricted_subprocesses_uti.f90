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

module restricted_subprocesses_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units, only: free_unit
  use format_defs, only: FMT_10, FMT_12
  use lorentz, only: vector4_t, vector3_moving, vector4_moving
  use particle_specifiers, only: new_prt_spec
  use process_libraries, only: process_library_t
  use resonances, only: resonance_info_t
  use resonances, only: resonance_history_t
  use resonances, only: resonance_history_set_t
  use state_matrices, only: FM_IGNORE_HELICITY
  use particles, only: particle_set_t
  use model_data, only: model_data_t
  use models, only: syntax_model_file_init, syntax_model_file_final
  use models, only: model_t
  use rng_base_ut, only: rng_test_factory_t
  use mci_base, only: mci_t
  use phs_base, only: phs_config_t
  use phs_forests, only: syntax_phs_forest_init, syntax_phs_forest_final
  use phs_wood, only: phs_wood_config_t
  use process_libraries, only: process_def_entry_t
  use process_libraries, only: process_component_def_t
  use prclib_stacks, only: prclib_entry_t
  use prc_core_def, only: prc_core_def_t
  use prc_omega, only: omega_def_t
  use process, only: process_t
  use instances, only: process_instance_t
  use process_stacks, only: process_entry_t
  use event_transforms, only: evt_trivial_t
  use resonance_insertion, only: evt_resonance_t
  use integrations, only: integrate_process
  use rt_data, only: rt_data_t

  use restricted_subprocesses

  implicit none
  private

  public :: restricted_subprocesses_1
  public :: restricted_subprocesses_2
  public :: restricted_subprocesses_3
  public :: restricted_subprocesses_4
  public :: restricted_subprocesses_5
  public :: restricted_subprocesses_6



  public :: prepare_resonance_test_library

contains

  subroutine restricted_subprocesses_1 (u)
    integer, intent(in) :: u
    type(rt_data_t) :: global
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(string_t) :: prc_name
    type(string_t), dimension(2) :: prt_in
    type(string_t), dimension(3) :: prt_out
    type(restricted_process_configuration_t) :: prc_config

    write (u, "(A)")  "* Test output: restricted_subprocesses_1"
    write (u, "(A)")  "*   Purpose: create subprocess list from resonances"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%select_model (var_str ("SM"))

    write (u, "(A)")  "* Create resonance history"
    write (u, "(A)")

    call res_info%init (3, -24, global%model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Create process configuration"
    write (u, "(A)")

    prc_name = "restricted_subprocesses_1_p"
    prt_in(1) = "e-"
    prt_in(2) = "e+"
    prt_out(1) = "d"
    prt_out(2) = "u"
    prt_out(3) = "W+"

    call prc_config%init_resonant_process (prc_name, &
         new_prt_spec (prt_in), new_prt_spec (prt_out), &
         res_history, global%model, global%var_list)

    call prc_config%write (u)

    write (u, *)
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: restricted_subprocesses_1"

  end subroutine restricted_subprocesses_1

  subroutine restricted_subprocesses_2 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global
    type(resonance_info_t) :: res_info
    type(resonance_history_t), dimension(2) :: res_history
    type(resonance_history_set_t) :: res_history_set
    type(string_t) :: libname
    type(string_t), dimension(2) :: prt_in
    type(string_t), dimension(3) :: prt_out
    type(resonant_subprocess_set_t) :: prc_set
    type(process_library_t), pointer :: lib
    logical :: exist

    write (u, "(A)")  "* Test output: restricted_subprocesses_2"
    write (u, "(A)")  "*   Purpose: create subprocess library from resonances"
    write (u, "(A)")

    call syntax_model_file_init ()

    call global%global_init ()
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%select_model (var_str ("SM"))

    write (u, "(A)")  "* Create resonance histories"
    write (u, "(A)")

    call res_info%init (3, -24, global%model, 5)
    call res_history(1)%add_resonance (res_info)
    call res_history(1)%write (u)

    call res_info%init (7, 23, global%model, 5)
    call res_history(2)%add_resonance (res_info)
    call res_history(2)%write (u)

    call res_history_set%init ()
    call res_history_set%enter (res_history(1))
    call res_history_set%enter (res_history(2))
    call res_history_set%freeze ()

    write (u, "(A)")
    write (u, "(A)")  "* Empty restricted subprocess set"
    write (u, "(A)")

    write (u, "(A,1x,L1)")  "active =", prc_set%is_active ()
    write (u, "(A)")

    call prc_set%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Fill restricted subprocess set"
    write (u, "(A)")

    libname = "restricted_subprocesses_2_p_R"
    prt_in(1) = "e-"
    prt_in(2) = "e+"
    prt_out(1) = "d"
    prt_out(2) = "u"
    prt_out(3) = "W+"

    call prc_set%init (1)
    call prc_set%fill_resonances (res_history_set, 1)
    call prc_set%create_library (libname, global, exist)
    if (.not. exist) then
       call prc_set%add_to_library (1, &
            new_prt_spec (prt_in), new_prt_spec (prt_out), &
            global)
    end if
    call prc_set%freeze_library (global)

    write (u, "(A,1x,L1)")  "active =", prc_set%is_active ()
    write (u, "(A)")

    call prc_set%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Queries"
    write (u, "(A)")

    write (u, "(A,1x,I0)")  "n_process =", prc_set%get_n_process ()
    write (u, "(A)")
    write (u, "(A,A,A)")  "libname = '", char (prc_set%get_libname ()), "'"
    write (u, "(A)")
    write (u, "(A,A,A)")  "proc_id(1) = '", char (prc_set%get_proc_id (1)), "'"
    write (u, "(A,A,A)")  "proc_id(2) = '", char (prc_set%get_proc_id (2)), "'"


    write (u, "(A)")
    write (u, "(A)")  "* Process library"
    write (u, "(A)")

    call prc_set%compile_library (global)

    lib => global%prclib_stack%get_library_ptr (libname)
    if (associated (lib))  call lib%write (u, libpath=.false.)

    write (u, *)
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: restricted_subprocesses_2"

  end subroutine restricted_subprocesses_2

  subroutine restricted_subprocesses_3 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global
    class(model_t), pointer :: model
    class(model_data_t), pointer :: model_data
    type(string_t) :: libname, libname_res
    type(string_t) :: procname
    type(process_component_def_t), pointer :: process_component_def
    type(prclib_entry_t), pointer :: lib_entry
    type(process_library_t), pointer :: lib
    logical :: exist
    type(process_t), pointer :: process
    type(process_instance_t), target :: process_instance
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(resonant_subprocess_set_t) :: prc_set
    type(particle_set_t) :: pset
    real(default) :: sqrts, mw, pp
    real(default), dimension(3) :: p3
    type(vector4_t), dimension(:), allocatable :: p
    real(default), dimension(:), allocatable :: m
    integer, dimension(:), allocatable :: pdg
    real(default), dimension(:), allocatable :: sqme
    logical, dimension(:), allocatable :: mask
    real(default) :: on_shell_limit
    integer, dimension(:), allocatable :: i_array
    real(default), dimension(:), allocatable :: prob_array
    type(evt_resonance_t), target :: evt_resonance
    integer :: i, u_dump

    write (u, "(A)")  "* Test output: restricted_subprocesses_3"
    write (u, "(A)")  "*   Purpose: handle process and resonance kinematics"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    call global%global_init ()
    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)
    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%set_log (var_str ("?resonance_history"), &
         .true., is_known = .true.)

    call global%select_model (var_str ("SM"))
    allocate (model)
    call model%init_instance (global%model)
    model_data => model

    libname = "restricted_subprocesses_3_lib"
    libname_res = "restricted_subprocesses_3_lib_res"
    procname = "restricted_subprocesses_3_p"

    write (u, "(A)")  "* Initialize process library and process"
    write (u, "(A)")

    allocate (lib_entry)
    call lib_entry%init (libname)
    lib => lib_entry%process_library_t
    call global%add_prclib (lib_entry)

    call prepare_resonance_test_library &
         (lib, libname, procname, model_data, global, u)

    call integrate_process (procname, global, &
         local_stack = .true., init_only = .true.)

    process => global%process_stack%get_process_ptr (procname)

    call process_instance%init (process)
    call process_instance%setup_event_data ()

    write (u, "(A)")
    write (u, "(A)")  "* Extract resonance history set"
    write (u, "(A)")

    call process%extract_resonance_history_set &
         (res_history_set(1), include_trivial=.true., i_component=1)
    call res_history_set(1)%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Build resonant-subprocess library"
    write (u, "(A)")

    call prc_set%init (1)
    call prc_set%fill_resonances (res_history_set(1), 1)

    process_component_def => process%get_component_def_ptr (1)
    call prc_set%create_library (libname_res, global, exist)
    if (.not. exist) then
       call prc_set%add_to_library (1, &
            process_component_def%get_prt_spec_in (), &
            process_component_def%get_prt_spec_out (), &
            global)
    end if
    call prc_set%freeze_library (global)
    call prc_set%compile_library (global)
    call prc_set%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Build particle set"
    write (u, "(A)")

    sqrts = global%get_rval (var_str ("sqrts"))
    mw = 80._default   ! deliberately slightly different from true mw
    pp = sqrt (sqrts**2 - 4 * mw**2) / 2

    allocate (pdg (5), p (5), m (5))
    pdg(1) = -11
    p(1) = vector4_moving (sqrts/2, sqrts/2, 3)
    m(1) = 0
    pdg(2) = 11
    p(2) = vector4_moving (sqrts/2,-sqrts/2, 3)
    m(2) = 0
    pdg(3) = 1
    p3(1) = pp/2
    p3(2) = mw/2
    p3(3) = 0
    p(3) = vector4_moving (sqrts/4, vector3_moving (p3))
    m(3) = 0
    p3(2) = -mw/2
    pdg(4) = -2
    p(4) = vector4_moving (sqrts/4, vector3_moving (p3))
    m(4) = 0
    pdg(5) = 24
    p(5) = vector4_moving (sqrts/2,-pp, 1)
    m(5) = mw

    call pset%init_direct (0, 2, 0, 0, 3, pdg, model)
    call pset%set_momentum (p, m**2)
    call pset%write (u, testflag=.true.)

     write (u, "(A)")
     write (u, "(A)")  "* Fill process instance"

    ! workflow from event_recalculate
    call process_instance%choose_mci (1)
    call process_instance%set_trace (pset, 1)
    call process_instance%recover &
         (1, 1, update_sqme=.true., recover_phs=.false.)
    call process_instance%evaluate_event_data (weight = 1._default)

    write (u, "(A)")
    write (u, "(A)")  "* Prepare resonant subprocesses"

    call prc_set%prepare_process_objects (global)
    call prc_set%prepare_process_instances (global)

    call evt_resonance%set_resonance_data (res_history_set)
    call evt_resonance%select_component (1)
    call prc_set%connect_transform (evt_resonance)
    call evt_resonance%connect (process_instance, model)

    call prc_set%fill_momenta ()

    write (u, "(A)")
    write (u, "(A)")  "* Show squared matrix element of master process,"
    write (u, "(A)")  "  should coincide with 2nd subprocess sqme"
    write (u, "(A)")

    write (u, "(1x,I0,1x," // FMT_12 // ")")  0, prc_set%get_master_sqme ()

    write (u, "(A)")
    write (u, "(A)")  "* Compute squared matrix elements &
         &of selected resonant subprocesses [1,2]"
    write (u, "(A)")

    call prc_set%evaluate_subprocess ([1,2])

    allocate (sqme (3), source = 0._default)
    call prc_set%get_subprocess_sqme (sqme)
    do i = 1, size (sqme)
       write (u, "(1x,I0,1x," // FMT_12 // ")")  i, sqme(i)
    end do
    deallocate (sqme)

    write (u, "(A)")
    write (u, "(A)")  "* Compute squared matrix elements &
         &of all resonant subprocesses"
    write (u, "(A)")

    call prc_set%evaluate_subprocess ([1,2,3])

    allocate (sqme (3), source = 0._default)
    call prc_set%get_subprocess_sqme (sqme)
    do i = 1, size (sqme)
       write (u, "(1x,I0,1x," // FMT_12 // ")")  i, sqme(i)
    end do
    deallocate (sqme)

    write (u, "(A)")
    write (u, "(A)")  "* Write process instances to file &
         &restricted_subprocesses_3_lib_res.dat"

    u_dump = free_unit ()
    open (unit = u_dump, file = "restricted_subprocesses_3_lib_res.dat", &
         action = "write", status = "replace")
    call prc_set%dump_instances (u_dump)
    close (u_dump)

    write (u, "(A)")
    write (u, "(A)")  "* Determine on-shell resonant subprocesses"
    write (u, "(A)")

    on_shell_limit = 0
    write (u, "(1x,A,1x," // FMT_10 // ")")  "on_shell_limit =", on_shell_limit
    call prc_set%set_on_shell_limit (on_shell_limit)
    call prc_set%determine_on_shell_histories (1, i_array)
    write (u, "(1x,A,9(1x,I0))")  "resonant =", i_array

    on_shell_limit = 0.1_default
    write (u, "(1x,A,1x," // FMT_10 // ")")  "on_shell_limit =", on_shell_limit
    call prc_set%set_on_shell_limit (on_shell_limit)
    call prc_set%determine_on_shell_histories (1, i_array)
    write (u, "(1x,A,9(1x,I0))")  "resonant =", i_array

    on_shell_limit = 10._default
    write (u, "(1x,A,1x," // FMT_10 // ")")  "on_shell_limit =", on_shell_limit
    call prc_set%set_on_shell_limit (on_shell_limit)
    call prc_set%determine_on_shell_histories (1, i_array)
    write (u, "(1x,A,9(1x,I0))")  "resonant =", i_array

    on_shell_limit = 10000._default
    write (u, "(1x,A,1x," // FMT_10 // ")")  "on_shell_limit =", on_shell_limit
    call prc_set%set_on_shell_limit (on_shell_limit)
    call prc_set%determine_on_shell_histories (1, i_array)
    write (u, "(1x,A,9(1x,I0))")  "resonant =", i_array

    write (u, "(A)")
    write (u, "(A)")  "* Compute probabilities for applicable resonances"
    write (u, "(A)")  "  and initialize the process selector"
    write (u, "(A)")  "  (The first number is the probability for background)"
    write (u, "(A)")

    on_shell_limit = 0
    write (u, "(1x,A,1x," // FMT_10 // ")")  "on_shell_limit =", on_shell_limit
    call prc_set%set_on_shell_limit (on_shell_limit)
    call prc_set%determine_on_shell_histories (1, i_array)
    call prc_set%compute_probabilities (prob_array)
    write (u, "(1x,A,9(1x,"// FMT_12 // "))")  "resonant =", prob_array
    call prc_set%write (u, testflag=.true.)
    write (u, *)

    on_shell_limit = 10._default
    write (u, "(1x,A,1x," // FMT_10 // ")")  "on_shell_limit =", on_shell_limit
    call prc_set%set_on_shell_limit (on_shell_limit)
    call prc_set%determine_on_shell_histories (1, i_array)
    call prc_set%compute_probabilities (prob_array)
    write (u, "(1x,A,9(1x,"// FMT_12 // "))")  "resonant =", prob_array
    call prc_set%write (u, testflag=.true.)
    write (u, *)

    on_shell_limit = 10000._default
    write (u, "(1x,A,1x," // FMT_10 // ")")  "on_shell_limit =", on_shell_limit
    call prc_set%set_on_shell_limit (on_shell_limit)
    call prc_set%determine_on_shell_histories (1, i_array)
    call prc_set%compute_probabilities (prob_array)
    write (u, "(1x,A,9(1x,"// FMT_12 // "))")  "resonant =", prob_array
    write (u, *)
    call prc_set%write (u, testflag=.true.)
    write (u, *)

    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_phs_forest_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: restricted_subprocesses_3"

  end subroutine restricted_subprocesses_3

  subroutine restricted_subprocesses_4 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global
    class(model_t), pointer :: model
    class(model_data_t), pointer :: model_data
    type(string_t) :: libname, libname_res
    type(string_t) :: procname
    type(process_component_def_t), pointer :: process_component_def
    type(prclib_entry_t), pointer :: lib_entry
    type(process_library_t), pointer :: lib
    logical :: exist
    type(process_t), pointer :: process
    type(process_instance_t), target :: process_instance
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(resonant_subprocess_set_t) :: prc_set
    type(particle_set_t) :: pset
    real(default) :: sqrts, mw, pp
    real(default), dimension(3) :: p3
    type(vector4_t), dimension(:), allocatable :: p
    real(default), dimension(:), allocatable :: m
    integer, dimension(:), allocatable :: pdg
    real(default) :: on_shell_limit
    type(evt_trivial_t), target :: evt_trivial
    type(evt_resonance_t), target :: evt_resonance
    real(default) :: probability
    integer :: i

    write (u, "(A)")  "* Test output: restricted_subprocesses_4"
    write (u, "(A)")  "*   Purpose: employ event transform"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    call global%global_init ()
    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)
    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%set_log (var_str ("?resonance_history"), &
         .true., is_known = .true.)

    call global%select_model (var_str ("SM"))
    allocate (model)
    call model%init_instance (global%model)
    model_data => model

    libname = "restricted_subprocesses_4_lib"
    libname_res = "restricted_subprocesses_4_lib_res"
    procname = "restricted_subprocesses_4_p"

    write (u, "(A)")  "* Initialize process library and process"
    write (u, "(A)")

    allocate (lib_entry)
    call lib_entry%init (libname)
    lib => lib_entry%process_library_t
    call global%add_prclib (lib_entry)

    call prepare_resonance_test_library &
         (lib, libname, procname, model_data, global, u)

    call integrate_process (procname, global, &
         local_stack = .true., init_only = .true.)

    process => global%process_stack%get_process_ptr (procname)

    call process_instance%init (process)
    call process_instance%setup_event_data ()

    write (u, "(A)")
    write (u, "(A)")  "* Extract resonance history set"

    call process%extract_resonance_history_set &
         (res_history_set(1), include_trivial=.false., i_component=1)

    write (u, "(A)")
    write (u, "(A)")  "* Build resonant-subprocess library"

    call prc_set%init (1)
    call prc_set%fill_resonances (res_history_set(1), 1)

    process_component_def => process%get_component_def_ptr (1)
    call prc_set%create_library (libname_res, global, exist)
    if (.not. exist) then
       call prc_set%add_to_library (1, &
            process_component_def%get_prt_spec_in (), &
            process_component_def%get_prt_spec_out (), &
            global)
    end if
    call prc_set%freeze_library (global)
    call prc_set%compile_library (global)

    write (u, "(A)")
    write (u, "(A)")  "* Build particle set"
    write (u, "(A)")

    sqrts = global%get_rval (var_str ("sqrts"))
    mw = 80._default   ! deliberately slightly different from true mw
    pp = sqrt (sqrts**2 - 4 * mw**2) / 2

    allocate (pdg (5), p (5), m (5))
    pdg(1) = -11
    p(1) = vector4_moving (sqrts/2, sqrts/2, 3)
    m(1) = 0
    pdg(2) = 11
    p(2) = vector4_moving (sqrts/2,-sqrts/2, 3)
    m(2) = 0
    pdg(3) = 1
    p3(1) = pp/2
    p3(2) = mw/2
    p3(3) = 0
    p(3) = vector4_moving (sqrts/4, vector3_moving (p3))
    m(3) = 0
    p3(2) = -mw/2
    pdg(4) = -2
    p(4) = vector4_moving (sqrts/4, vector3_moving (p3))
    m(4) = 0
    pdg(5) = 24
    p(5) = vector4_moving (sqrts/2,-pp, 1)
    m(5) = mw

    call pset%init_direct (0, 2, 0, 0, 3, pdg, model)
    call pset%set_momentum (p, m**2)

    write (u, "(A)")  "* Fill process instance"
    write (u, "(A)")

    ! workflow from event_recalculate
    call process_instance%choose_mci (1)
    call process_instance%set_trace (pset, 1)
    call process_instance%recover &
         (1, 1, update_sqme=.true., recover_phs=.false.)
    call process_instance%evaluate_event_data (weight = 1._default)

    write (u, "(A)")  "* Prepare resonant subprocesses"
    write (u, "(A)")

    call prc_set%prepare_process_objects (global)
    call prc_set%prepare_process_instances (global)

    write (u, "(A)")  "* Fill trivial event transform (deliberately w/o color)"
    write (u, "(A)")

    call evt_trivial%connect (process_instance, model)
    call evt_trivial%set_particle_set (pset, 1, 1)
    call evt_trivial%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize resonance-insertion event transform"
    write (u, "(A)")

    evt_trivial%next => evt_resonance
    evt_resonance%previous => evt_trivial

    call evt_resonance%set_resonance_data (res_history_set)
    call evt_resonance%select_component (1)
    call evt_resonance%connect (process_instance, model)
    call prc_set%connect_transform (evt_resonance)
    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute probabilities for applicable resonances"
    write (u, "(A)")  "  and initialize the process selector"
    write (u, "(A)")

    on_shell_limit = 10._default
    write (u, "(1x,A,1x," // FMT_10 // ")")  "on_shell_limit =", on_shell_limit
    call evt_resonance%set_on_shell_limit (on_shell_limit)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate resonance-insertion event transform"
    write (u, "(A)")

    call evt_resonance%prepare_new_event (1, 1)
    call evt_resonance%generate_weighted (probability)
    call evt_resonance%make_particle_set (1, .false.)

    call evt_resonance%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_phs_forest_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: restricted_subprocesses_4"

  end subroutine restricted_subprocesses_4

  subroutine restricted_subprocesses_5 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global
    class(model_t), pointer :: model
    class(model_data_t), pointer :: model_data
    type(string_t) :: libname, libname_res
    type(string_t) :: procname
    type(process_component_def_t), pointer :: process_component_def
    type(prclib_entry_t), pointer :: lib_entry
    type(process_library_t), pointer :: lib
    logical :: exist
    type(process_t), pointer :: process
    type(process_instance_t), target :: process_instance
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(resonant_subprocess_set_t) :: prc_set
    type(particle_set_t) :: pset
    real(default) :: sqrts, mw, pp
    real(default), dimension(3) :: p3
    type(vector4_t), dimension(:), allocatable :: p
    real(default), dimension(:), allocatable :: m
    integer, dimension(:), allocatable :: pdg
    real(default) :: on_shell_limit
    real(default) :: on_shell_turnoff
    type(evt_trivial_t), target :: evt_trivial
    type(evt_resonance_t), target :: evt_resonance
    real(default) :: probability
    integer :: i

    write (u, "(A)")  "* Test output: restricted_subprocesses_5"
    write (u, "(A)")  "*   Purpose: employ event transform &
         &with gaussian turnoff"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    call global%global_init ()
    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)
    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%set_log (var_str ("?resonance_history"), &
         .true., is_known = .true.)

    call global%select_model (var_str ("SM"))
    allocate (model)
    call model%init_instance (global%model)
    model_data => model

    libname = "restricted_subprocesses_5_lib"
    libname_res = "restricted_subprocesses_5_lib_res"
    procname = "restricted_subprocesses_5_p"

    write (u, "(A)")  "* Initialize process library and process"
    write (u, "(A)")

    allocate (lib_entry)
    call lib_entry%init (libname)
    lib => lib_entry%process_library_t
    call global%add_prclib (lib_entry)

    call prepare_resonance_test_library &
         (lib, libname, procname, model_data, global, u)

    call integrate_process (procname, global, &
         local_stack = .true., init_only = .true.)

    process => global%process_stack%get_process_ptr (procname)

    call process_instance%init (process)
    call process_instance%setup_event_data ()

    write (u, "(A)")
    write (u, "(A)")  "* Extract resonance history set"

    call process%extract_resonance_history_set &
         (res_history_set(1), include_trivial=.false., i_component=1)

    write (u, "(A)")
    write (u, "(A)")  "* Build resonant-subprocess library"

    call prc_set%init (1)
    call prc_set%fill_resonances (res_history_set(1), 1)

    process_component_def => process%get_component_def_ptr (1)
    call prc_set%create_library (libname_res, global, exist)
    if (.not. exist) then
       call prc_set%add_to_library (1, &
            process_component_def%get_prt_spec_in (), &
            process_component_def%get_prt_spec_out (), &
            global)
    end if
    call prc_set%freeze_library (global)
    call prc_set%compile_library (global)

    write (u, "(A)")
    write (u, "(A)")  "* Build particle set"
    write (u, "(A)")

    sqrts = global%get_rval (var_str ("sqrts"))
    mw = 80._default   ! deliberately slightly different from true mw
    pp = sqrt (sqrts**2 - 4 * mw**2) / 2

    allocate (pdg (5), p (5), m (5))
    pdg(1) = -11
    p(1) = vector4_moving (sqrts/2, sqrts/2, 3)
    m(1) = 0
    pdg(2) = 11
    p(2) = vector4_moving (sqrts/2,-sqrts/2, 3)
    m(2) = 0
    pdg(3) = 1
    p3(1) = pp/2
    p3(2) = mw/2
    p3(3) = 0
    p(3) = vector4_moving (sqrts/4, vector3_moving (p3))
    m(3) = 0
    p3(2) = -mw/2
    pdg(4) = -2
    p(4) = vector4_moving (sqrts/4, vector3_moving (p3))
    m(4) = 0
    pdg(5) = 24
    p(5) = vector4_moving (sqrts/2,-pp, 1)
    m(5) = mw

    call pset%init_direct (0, 2, 0, 0, 3, pdg, model)
    call pset%set_momentum (p, m**2)

    write (u, "(A)")  "* Fill process instance"
    write (u, "(A)")

    ! workflow from event_recalculate
    call process_instance%choose_mci (1)
    call process_instance%set_trace (pset, 1)
    call process_instance%recover &
         (1, 1, update_sqme=.true., recover_phs=.false.)
    call process_instance%evaluate_event_data (weight = 1._default)

    write (u, "(A)")  "* Prepare resonant subprocesses"
    write (u, "(A)")

    call prc_set%prepare_process_objects (global)
    call prc_set%prepare_process_instances (global)

    write (u, "(A)")  "* Fill trivial event transform (deliberately w/o color)"
    write (u, "(A)")

    call evt_trivial%connect (process_instance, model)
    call evt_trivial%set_particle_set (pset, 1, 1)
    call evt_trivial%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize resonance-insertion event transform"
    write (u, "(A)")

    evt_trivial%next => evt_resonance
    evt_resonance%previous => evt_trivial

    call evt_resonance%set_resonance_data (res_history_set)
    call evt_resonance%select_component (1)
    call evt_resonance%connect (process_instance, model)
    call prc_set%connect_transform (evt_resonance)
    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute probabilities for applicable resonances"
    write (u, "(A)")  "  and initialize the process selector"
    write (u, "(A)")

    on_shell_limit = 10._default
    write (u, "(1x,A,1x," // FMT_10 // ")")  "on_shell_limit   =", &
         on_shell_limit
    call evt_resonance%set_on_shell_limit (on_shell_limit)

    on_shell_turnoff = 1._default
    write (u, "(1x,A,1x," // FMT_10 // ")")  "on_shell_turnoff =", &
         on_shell_turnoff
    call evt_resonance%set_on_shell_turnoff (on_shell_turnoff)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate resonance-insertion event transform"
    write (u, "(A)")

    call evt_resonance%prepare_new_event (1, 1)
    call evt_resonance%generate_weighted (probability)
    call evt_resonance%make_particle_set (1, .false.)

    call evt_resonance%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_phs_forest_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: restricted_subprocesses_5"

  end subroutine restricted_subprocesses_5

  subroutine restricted_subprocesses_6 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global
    class(model_t), pointer :: model
    class(model_data_t), pointer :: model_data
    type(string_t) :: libname, libname_res
    type(string_t) :: procname
    type(process_component_def_t), pointer :: process_component_def
    type(prclib_entry_t), pointer :: lib_entry
    type(process_library_t), pointer :: lib
    logical :: exist
    type(process_t), pointer :: process
    type(process_instance_t), target :: process_instance
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(resonant_subprocess_set_t) :: prc_set
    type(particle_set_t) :: pset
    real(default) :: sqrts, mw, pp
    real(default), dimension(3) :: p3
    type(vector4_t), dimension(:), allocatable :: p
    real(default), dimension(:), allocatable :: m
    integer, dimension(:), allocatable :: pdg
    real(default) :: on_shell_limit
    real(default) :: background_factor
    type(evt_trivial_t), target :: evt_trivial
    type(evt_resonance_t), target :: evt_resonance
    real(default) :: probability
    integer :: i

    write (u, "(A)")  "* Test output: restricted_subprocesses_6"
    write (u, "(A)")  "*   Purpose: employ event transform &
         &with background switched off"
    write (u, "(A)")

    call syntax_model_file_init ()
    call syntax_phs_forest_init ()

    call global%global_init ()
    call global%append_log (&
         var_str ("?rebuild_phase_space"), .true., intrinsic = .true.)
    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)
    call global%set_int (var_str ("seed"), &
         0, is_known = .true.)
    call global%set_real (var_str ("sqrts"),&
         1000._default, is_known = .true.)
    call global%set_log (var_str ("?resonance_history"), &
         .true., is_known = .true.)

    call global%select_model (var_str ("SM"))
    allocate (model)
    call model%init_instance (global%model)
    model_data => model

    libname = "restricted_subprocesses_6_lib"
    libname_res = "restricted_subprocesses_6_lib_res"
    procname = "restricted_subprocesses_6_p"

    write (u, "(A)")  "* Initialize process library and process"
    write (u, "(A)")

    allocate (lib_entry)
    call lib_entry%init (libname)
    lib => lib_entry%process_library_t
    call global%add_prclib (lib_entry)

    call prepare_resonance_test_library &
         (lib, libname, procname, model_data, global, u)

    call integrate_process (procname, global, &
         local_stack = .true., init_only = .true.)

    process => global%process_stack%get_process_ptr (procname)

    call process_instance%init (process)
    call process_instance%setup_event_data ()

    write (u, "(A)")
    write (u, "(A)")  "* Extract resonance history set"

    call process%extract_resonance_history_set &
         (res_history_set(1), include_trivial=.false., i_component=1)

    write (u, "(A)")
    write (u, "(A)")  "* Build resonant-subprocess library"

    call prc_set%init (1)
    call prc_set%fill_resonances (res_history_set(1), 1)

    process_component_def => process%get_component_def_ptr (1)
    call prc_set%create_library (libname_res, global, exist)
    if (.not. exist) then
       call prc_set%add_to_library (1, &
            process_component_def%get_prt_spec_in (), &
            process_component_def%get_prt_spec_out (), &
            global)
    end if
    call prc_set%freeze_library (global)
    call prc_set%compile_library (global)

    write (u, "(A)")
    write (u, "(A)")  "* Build particle set"
    write (u, "(A)")

    sqrts = global%get_rval (var_str ("sqrts"))
    mw = 80._default   ! deliberately slightly different from true mw
    pp = sqrt (sqrts**2 - 4 * mw**2) / 2

    allocate (pdg (5), p (5), m (5))
    pdg(1) = -11
    p(1) = vector4_moving (sqrts/2, sqrts/2, 3)
    m(1) = 0
    pdg(2) = 11
    p(2) = vector4_moving (sqrts/2,-sqrts/2, 3)
    m(2) = 0
    pdg(3) = 1
    p3(1) = pp/2
    p3(2) = mw/2
    p3(3) = 0
    p(3) = vector4_moving (sqrts/4, vector3_moving (p3))
    m(3) = 0
    p3(2) = -mw/2
    pdg(4) = -2
    p(4) = vector4_moving (sqrts/4, vector3_moving (p3))
    m(4) = 0
    pdg(5) = 24
    p(5) = vector4_moving (sqrts/2,-pp, 1)
    m(5) = mw

    call pset%init_direct (0, 2, 0, 0, 3, pdg, model)
    call pset%set_momentum (p, m**2)

    write (u, "(A)")  "* Fill process instance"
    write (u, "(A)")

    ! workflow from event_recalculate
    call process_instance%choose_mci (1)
    call process_instance%set_trace (pset, 1)
    call process_instance%recover &
         (1, 1, update_sqme=.true., recover_phs=.false.)
    call process_instance%evaluate_event_data (weight = 1._default)

    write (u, "(A)")  "* Prepare resonant subprocesses"
    write (u, "(A)")

    call prc_set%prepare_process_objects (global)
    call prc_set%prepare_process_instances (global)

    write (u, "(A)")  "* Fill trivial event transform (deliberately w/o color)"
    write (u, "(A)")

    call evt_trivial%connect (process_instance, model)
    call evt_trivial%set_particle_set (pset, 1, 1)
    call evt_trivial%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize resonance-insertion event transform"
    write (u, "(A)")

    evt_trivial%next => evt_resonance
    evt_resonance%previous => evt_trivial

    call evt_resonance%set_resonance_data (res_history_set)
    call evt_resonance%select_component (1)
    call evt_resonance%connect (process_instance, model)
    call prc_set%connect_transform (evt_resonance)
    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute probabilities for applicable resonances"
    write (u, "(A)")  "  and initialize the process selector"
    write (u, "(A)")

    on_shell_limit = 10._default
    write (u, "(1x,A,1x," // FMT_10 // ")") &
         "on_shell_limit    =", on_shell_limit
    call evt_resonance%set_on_shell_limit (on_shell_limit)

    background_factor = 0
    write (u, "(1x,A,1x," // FMT_10 // ")") &
         "background_factor =", background_factor
    call evt_resonance%set_background_factor (background_factor)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate resonance-insertion event transform"
    write (u, "(A)")

    call evt_resonance%prepare_new_event (1, 1)
    call evt_resonance%generate_weighted (probability)
    call evt_resonance%make_particle_set (1, .false.)

    call evt_resonance%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_phs_forest_final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: restricted_subprocesses_6"

  end subroutine restricted_subprocesses_6


  subroutine prepare_resonance_test_library &
       (lib, libname, procname, model, global, u)
    type(process_library_t), target, intent(out) :: lib
    type(string_t), intent(in) :: libname
    type(string_t), intent(in) :: procname
    class(model_data_t), intent(in), pointer :: model
    type(rt_data_t), intent(in), target :: global
    integer, intent(in) :: u
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    class(prc_core_def_t), allocatable :: def
    type(process_def_entry_t), pointer :: entry

    call lib%init (libname)

    allocate (prt_in (2), prt_out (3))
    prt_in = [var_str ("e+"), var_str ("e-")]
    prt_out = [var_str ("d"), var_str ("ubar"), var_str ("W+")]

    allocate (omega_def_t :: def)
    select type (def)
    type is (omega_def_t)
       call def%init (model%get_name (), prt_in, prt_out, &
            ovm=.false., ufo=.false.)
    end select
    allocate (entry)
    call entry%init (procname, &
         model_name = model%get_name (), &
         n_in = 2, n_components = 1, &
         requires_resonances = .true.)
    call entry%import_component (1, n_out = size (prt_out), &
         prt_in  = new_prt_spec (prt_in), &
         prt_out = new_prt_spec (prt_out), &
         method  = var_str ("omega"), &
         variant = def)
    call entry%write (u)

    call lib%append (entry)

    call lib%configure (global%os_data)
    call lib%write_makefile (global%os_data, force = .true., verbose = .false.)
    call lib%clean (global%os_data, distclean = .false.)
    call lib%write_driver (force = .true.)
    call lib%load (global%os_data)

  end subroutine prepare_resonance_test_library


end module restricted_subprocesses_uti

