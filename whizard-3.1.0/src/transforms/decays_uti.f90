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

module decays_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use sm_qcd
  use model_data
  use models
  use state_matrices, only: FM_IGNORE_HELICITY
  use interactions, only: reset_interaction_counter
  use flavors
  use process_libraries
  use rng_base
  use mci_base
  use mci_midpoint
  use phs_base
  use phs_single
  use prc_core
  use prc_test, only: prc_test_create_library
  use process, only: process_t
  use instances, only: process_instance_t
  use process_stacks

  use decays

  use rng_base_ut, only: rng_test_t, rng_test_factory_t

  implicit none
  private

  public :: prepare_testbed

  public :: decays_1
  public :: decays_2
  public :: decays_3
  public :: decays_4
  public :: decays_5
  public :: decays_6

contains

  subroutine decays_1 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(model_data_t), target :: model
    type(flavor_t) :: flv_h
    type(flavor_t), dimension(2,1) :: flv_hbb, flv_hgg
    type(unstable_config_t), allocatable :: unstable

    write (u, "(A)")  "* Test output: decays_1"
    write (u, "(A)")  "*   Purpose: Set up branching and decay configuration"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment"
    write (u, "(A)")

    call os_data%init ()
    call model%init_sm_test ()

    call flv_h%init (25, model)
    call flv_hbb(:,1)%init ([5, -5], model)
    call flv_hgg(:,1)%init ([22, 22], model)

    write (u, "(A)")  "* Set up branching and decay"
    write (u, "(A)")

    allocate (unstable)
    call unstable%init (flv_h)
    call unstable%init_decays ([var_str ("h_bb"), var_str ("h_gg")], model)

    call unstable%init_test_case1 &
         (1, flv_hbb, 1.234e-3_default, .02_default, model)

    call unstable%init_test_case1 &
         (2, flv_hgg, 3.085e-4_default, .08_default, model)

    call unstable%compute ()
    call unstable%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call unstable%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: decays_1"

  end subroutine decays_1

  subroutine decays_2 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(model_data_t), target :: model
    type(flavor_t) :: flv_h, flv_wp, flv_wm
    type(flavor_t), dimension(2,1) :: flv_hww, flv_wud, flv_wen
    type(unstable_config_t), allocatable :: unstable

    write (u, "(A)")  "* Test output: decays_2"
    write (u, "(A)")  "*   Purpose: Set up cascade branching"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment"
    write (u, "(A)")

    call os_data%init ()
    call model%init_sm_test ()

    call model%set_unstable (25, [var_str ("h_ww")])
    call model%set_unstable (24, [var_str ("w_ud"), var_str ("w_en")])

    call flv_h%init (25, model)
    call flv_hww(:,1)%init ([24, -24], model)
    call flv_wp%init (24, model)
    call flv_wm%init (-24, model)
    call flv_wud(:,1)%init ([2, -1], model)
    call flv_wen(:,1)%init ([-11, 12], model)


    write (u, "(A)")  "* Set up branching and decay"
    write (u, "(A)")

    allocate (unstable)
    call unstable%init (flv_h, set_decays=.true., model=model)

    call unstable%init_test_case2 (flv_hww, flv_wud, flv_wen, model)

    call unstable%compute ()
    call unstable%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call unstable%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: decays_2"

  end subroutine decays_2

  subroutine decays_3 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    type(process_library_t), target :: lib
    type(string_t) :: prefix
    type(string_t) :: procname2
    type(process_stack_t) :: process_stack
    type(process_t), pointer :: process
    type(unstable_config_t), allocatable :: unstable
    type(flavor_t) :: flv

    write (u, "(A)")  "* Test output: decays_3"
    write (u, "(A)")  "*   Purpose: Connect a decay configuration &
         &with a process"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment and integrate process"
    write (u, "(A)")

    call os_data%init ()

    prefix = "decays_3"
    call prepare_testbed &
         (lib, process_stack, prefix, os_data, &
         scattering=.false., decay=.true., decay_rest_frame=.false.)

    procname2 = prefix // "_d"
    process => process_stack%get_process_ptr (procname2)
    model => process%get_model_ptr ()
    call process%write (.false., u)

    write (u, "(A)")
    write (u, "(A)")  "* Set up branching and decay"
    write (u, "(A)")

    call flv%init (25, model)

    allocate (unstable)
    call unstable%init (flv)
    call unstable%init_decays ([procname2], model)

    write (u, "(A)")  "* Connect decay with process object"
    write (u, "(A)")

    call unstable%connect_decay (1, process, model)

    call unstable%compute ()
    call unstable%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call unstable%final ()
    call process_stack%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: decays_3"

  end subroutine decays_3

  subroutine decays_4 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    type(process_library_t), target :: lib
    type(string_t) :: prefix, procname2
    class(rng_t), allocatable :: rng
    type(process_stack_t) :: process_stack
    type(process_t), pointer :: process
    type(unstable_config_t), allocatable, target :: unstable
    type(flavor_t) :: flv
    type(unstable_t), allocatable :: instance

    write (u, "(A)")  "* Test output: decays_4"
    write (u, "(A)")  "*   Purpose: Create a decay process and evaluate &
         &an instance"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment, process, &
         &and decay configuration"
    write (u, "(A)")

    call os_data%init ()

    prefix = "decays_4"
    call prepare_testbed &
         (lib, process_stack, prefix, os_data, &
         scattering=.false., decay=.true., decay_rest_frame = .false.)

    procname2 = prefix // "_d"
    process => process_stack%get_process_ptr (procname2)
    model => process%get_model_ptr ()

    call flv%init (25, model)

    allocate (unstable)
    call unstable%init (flv)
    call unstable%init_decays ([procname2], model)

    call model%set_unstable (25, [procname2])

    call unstable%connect_decay (1, process, model)

    call unstable%compute ()

    allocate (rng_test_t :: rng)

    allocate (instance)
    call instance%init (unstable)
    call instance%import_rng (rng)

    call instance%select_chain ()
    call instance%generate ()
    call instance%write (u)

    write (u, *)
    call instance%write_process_instances (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call instance%final ()
    call process_stack%final ()
    call unstable%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: decays_4"

  end subroutine decays_4

  subroutine decays_5 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    type(process_library_t), target :: lib
    type(string_t) :: prefix, procname1, procname2
    type(process_stack_t) :: process_stack
    type(process_t), pointer :: process
    type(process_instance_t), allocatable, target :: process_instance
    type(decay_root_config_t), target :: decay_root_config
    type(decay_root_t) :: decay_root
    type(decay_chain_t) :: decay_chain

    write (u, "(A)")  "* Test output: decays_5"
    write (u, "(A)")  "*   Purpose: Handle a process with subsequent decays"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment and parent process"
    write (u, "(A)")

    call os_data%init ()

    prefix = "decays_5"
    procname1 = prefix // "_p"
    procname2 = prefix // "_d"
    call prepare_testbed &
         (lib, process_stack, prefix, os_data, &
         scattering=.true., decay=.true.)

    write (u, "(A)")  "* Initialize decay process"
    write (u, "(A)")

    process => process_stack%get_process_ptr (procname1)
    model => process%get_model_ptr ()
    call model%set_unstable (25, [procname2])

    write (u, "(A)")  "* Initialize decay tree configuration"
    write (u, "(A)")

    call decay_root_config%connect (process, model, process_stack)
    call decay_root_config%compute ()
    call decay_root_config%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize decay tree"

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data ()
    call process_instance%init_simulation (1)

    call decay_root%init (decay_root_config, process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Select decay chain"
    write (u, "(A)")

    call decay_root%set_mci (1)
    !!! Not yet implemented; there is only one term anyway:
    ! call process_instance%select_i_term (decay_root%selected_term)
    call decay_root%set_term (1)
    call decay_root%select_chain ()

    call decay_chain%build (decay_root)

    call decay_root%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate event"
    write (u, "(A)")

    call process_instance%generate_unweighted_event (decay_root%get_mci ())
    call process_instance%evaluate_event_data ()

    call decay_root%generate ()

    call pacify (decay_root)

    write (u, "(A)")  "* Process instances"
    write (u, "(A)")

    call decay_root%write_process_instances (u)

    write (u, "(A)")
    write (u, "(A)")  "* Generate decay chain"
    write (u, "(A)")

    call decay_chain%evaluate ()
    call decay_chain%write (u)

    write (u, *)
    write (u, "(A,ES19.12)")  "chain probability =", &
         decay_chain%get_probability ()

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call decay_chain%final ()
    call decay_root%final ()
    call decay_root_config%final ()
    call process_instance%final ()
    deallocate (process_instance)

    call process_stack%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: decays_5"

  end subroutine decays_5

  subroutine decays_6 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    class(model_data_t), pointer :: model
    type(process_library_t), target :: lib
    type(string_t) :: prefix, procname1, procname2
    type(process_stack_t) :: process_stack
    type(process_t), pointer :: process
    type(process_instance_t), allocatable, target :: process_instance
    type(evt_decay_t), target :: evt_decay
    integer :: factorization_mode
    logical :: keep_correlations

    write (u, "(A)")  "* Test output: decays_6"
    write (u, "(A)")  "*   Purpose: Handle a process with subsequent decays"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment and parent process"
    write (u, "(A)")

    call os_data%init ()

    prefix = "decays_6"
    procname1 = prefix // "_p"
    procname2 = prefix // "_d"
    call prepare_testbed &
         (lib, process_stack, prefix, os_data, &
         scattering=.true., decay=.true.)

    write (u, "(A)")  "* Initialize decay process"

    process => process_stack%get_process_ptr (procname1)
    model => process%get_model_ptr ()
    call model%set_unstable (25, [procname2])

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data ()
    call process_instance%init_simulation (1)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize decay object"

    call evt_decay%connect (process_instance, model, process_stack)

    write (u, "(A)")
    write (u, "(A)")  "* Generate scattering event"

    call process_instance%generate_unweighted_event (1)
    call process_instance%evaluate_event_data ()

    write (u, "(A)")
    write (u, "(A)")  "* Select decay chain and generate event"
    write (u, "(A)")

    call evt_decay%prepare_new_event (1, 1)
    call evt_decay%generate_unweighted ()

    factorization_mode = FM_IGNORE_HELICITY
    keep_correlations = .false.
    call evt_decay%make_particle_set (factorization_mode, keep_correlations)

    call evt_decay%write (u, verbose = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_decay%final ()
    call process_instance%final ()
    deallocate (process_instance)

    call process_stack%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: decays_6"

  end subroutine decays_6


  subroutine prepare_testbed &
       (lib, process_stack, prefix, os_data, &
        scattering, decay, decay_rest_frame)
    type(process_library_t), intent(out), target :: lib
    type(process_stack_t), intent(out) :: process_stack
    type(string_t), intent(in) :: prefix
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: scattering, decay
    logical, intent(in), optional :: decay_rest_frame

    type(model_t), target :: model
    type(model_t), target :: model_copy
    type(string_t) :: libname, procname1, procname2
    type(process_entry_t), pointer :: process
    type(process_instance_t), allocatable, target :: process_instance
    class(phs_config_t), allocatable :: phs_config_template
    type(field_data_t), pointer :: field_data
    real(default) :: sqrts

    libname = prefix // "_lib"
    procname1 = prefix // "_p"
    procname2 = prefix // "_d"

    call model%init_test ()
    call model%set_par (var_str ("ff"), 0.4_default)
    call model%set_par (var_str ("mf"), &
         model%get_real (var_str ("ff")) * model%get_real (var_str ("ms")))

    if (scattering .and. decay) then
       field_data => model%get_field_ptr (25)
       call field_data%set (p_is_stable = .false.)
    end if

    call prc_test_create_library (libname, lib, &
         scattering = .true., decay = .true., &
         procname1 = procname1, procname2 = procname2)

    call reset_interaction_counter ()

    allocate (phs_single_config_t :: phs_config_template)

    if (scattering) then

       call model_copy%init (model%get_name (), &
            model%get_n_real (), &
            model%get_n_complex (), &
            model%get_n_field (), &
            model%get_n_vtx ())
       call model_copy%copy_from (model)

       allocate (process)
       call process%init (procname1, lib, os_data, model_copy)
       call process%setup_test_cores ()
       call process%init_components (phs_config_template)
       sqrts = 1000
       call process%setup_beams_sqrts (sqrts, i_core = 1)
       call process%configure_phs ()
       call process%setup_mci (dispatch_mci_test_midpoint)
       call process%setup_terms ()

       allocate (process_instance)
       call process_instance%init (process%process_t)
       call process_instance%integrate (1, n_it = 1, n_calls = 100)
       call process%final_integration (1)
       call process_instance%final ()
       deallocate (process_instance)

       call process%prepare_simulation (1)
       call process_stack%push (process)
    end if

    if (decay) then
       call model_copy%init (model%get_name (), &
            model%get_n_real (), &
            model%get_n_complex (), &
            model%get_n_field (), &
            model%get_n_vtx ())
       call model_copy%copy_from (model)

       allocate (process)
       call process%init (procname2, lib, os_data, model_copy)
       call process%setup_test_cores ()
       call process%init_components (phs_config_template)
       if (present (decay_rest_frame)) then
          call process%setup_beams_decay (rest_frame = decay_rest_frame, i_core = 1)
       else
          call process%setup_beams_decay (rest_frame = .not. scattering, i_core = 1)
       end if
       call process%configure_phs ()
       call process%setup_mci (dispatch_mci_test_midpoint)
       call process%setup_terms ()

       allocate (process_instance)
       call process_instance%init (process%process_t)
       call process_instance%integrate (1, n_it=1, n_calls=100)
       call process%final_integration (1)
       call process_instance%final ()
       deallocate (process_instance)

       call process%prepare_simulation (1)
       call process_stack%push (process)
    end if

    call model%final ()
    call model_copy%final ()

  end subroutine prepare_testbed

  subroutine dispatch_mci_test_midpoint (mci, var_list, process_id, is_nlo)
    use variables, only: var_list_t
    class(mci_t), allocatable, intent(out) :: mci
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_id
    logical, intent(in), optional :: is_nlo
    allocate (mci_midpoint_t :: mci)
  end subroutine dispatch_mci_test_midpoint


end module decays_uti

