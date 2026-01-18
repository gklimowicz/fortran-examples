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

submodule (restricted_subprocesses) restricted_subprocesses_s

  use diagnostics, only: msg_message, msg_fatal, msg_bug
  use diagnostics, only: signal_is_pending
  use io_units, only: given_output_unit
  use format_defs, only: FMT_14, FMT_19
  use string_utils, only: str
  use process_libraries, only: process_component_def_t
  use process_libraries, only: process_library_t
  use process_libraries, only: STAT_ACTIVE
  use prclib_stacks, only: prclib_entry_t
  use compilations, only: compile_library
  use integrations, only: integrate_process

  implicit none

contains

  module subroutine init_resonant_process &
       (prc_config, prc_name, prt_in, prt_out, res_history, model, var_list)
    class(restricted_process_configuration_t), intent(out) :: prc_config
    type(string_t), intent(in) :: prc_name
    type(prt_spec_t), dimension(:), intent(in) :: prt_in
    type(prt_spec_t), dimension(:), intent(in) :: prt_out
    type(resonance_history_t), intent(in) :: res_history
    type(model_t), intent(in), target :: model
    type(var_list_t), intent(in), target :: var_list
    type(model_t), pointer :: local_model
    type(var_list_t) :: local_var_list
    allocate (local_model)
    call local_model%init_instance (model)
    call local_var_list%link (var_list)
    call local_var_list%append_string (var_str ("$model_name"), &
         sval = local_model%get_name (), &
         intrinsic=.true.)
    call local_var_list%append_string (var_str ("$method"), &
         sval = var_str ("omega"), &
         intrinsic=.true.)
    call local_var_list%append_string (var_str ("$restrictions"), &
         sval = res_history%as_omega_string (size (prt_in)), &
         intrinsic = .true.)
    call local_var_list%append_log (var_str ("?resonance_history"), &
         lval = .false., &
         intrinsic = .true.)
    call prc_config%init (prc_name, size (prt_in), 1, &
         local_model, local_var_list)
    call prc_config%setup_component (1, &
         prt_in, prt_out, &
         local_model, local_var_list)
  end subroutine init_resonant_process

  module subroutine resonant_subprocess_set_write (prc_set, unit, testflag)
    class(resonant_subprocess_set_t), intent(in) :: prc_set
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    logical :: truncate
    integer :: u, i
    u = given_output_unit (unit)
    truncate = .false.;  if (present (testflag))  truncate = testflag
    write (u, "(1x,A)")  "Resonant subprocess set:"
    if (allocated (prc_set%n_history)) then
       if (any (prc_set%n_history > 0)) then
          do i = 1, size (prc_set%n_history)
             if (prc_set%n_history(i) > 0) then
                write (u, "(1x,A,I0)")  "Component #", i
                call prc_set%res_history_set(i)%write (u, indent=1)
             end if
          end do
          if (prc_set%lib_active) then
             write (u, "(3x,A,A,A)")  "Process library = '", &
                  char (prc_set%libname), "'"
          else
             write (u, "(3x,A)")  "Process library: [inactive]"
          end if
          if (associated (prc_set%evt)) then
             if (truncate) then
                write (u, "(3x,A,1x," // FMT_14 // ")") &
                     "Process sqme =", prc_set%get_master_sqme ()
             else
                write (u, "(3x,A,1x," // FMT_19 // ")") &
                     "Process sqme =", prc_set%get_master_sqme ()
             end if
          end if
          if (associated (prc_set%evt)) then
             write (u, "(3x,A)")  "Event transform: associated"
             write (u, "(2x)", advance="no")
             call prc_set%evt%write_selector (u, testflag)
          else
             write (u, "(3x,A)")  "Event transform: not associated"
          end if
       else
          write (u, "(2x,A)")  "[empty]"
       end if
    else
       write (u, "(3x,A)")  "[not allocated]"
    end if
  end subroutine resonant_subprocess_set_write

  module subroutine resonant_subprocess_set_init (prc_set, n_component)
    class(resonant_subprocess_set_t), intent(out) :: prc_set
    integer, intent(in) :: n_component
    allocate (prc_set%res_history_set (n_component))
    allocate (prc_set%n_history (n_component), source = 0)
  end subroutine resonant_subprocess_set_init

  module subroutine resonant_subprocess_set_fill_resonances (prc_set, &
       res_history_set, i_component)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    type(resonance_history_set_t), intent(in) :: res_history_set
    integer, intent(in) :: i_component
    prc_set%n_history(i_component) = res_history_set%get_n_history ()
    if (prc_set%n_history(i_component) > 0) then
       prc_set%res_history_set(i_component) = res_history_set
    else
       call prc_set%res_history_set(i_component)%init (initial_size = 0)
       call prc_set%res_history_set(i_component)%freeze ()
    end if
  end subroutine resonant_subprocess_set_fill_resonances

  module function resonant_subprocess_set_get_resonance_history_set &
       (prc_set) result (res_history_set)
    class(resonant_subprocess_set_t), intent(in) :: prc_set
    type(resonance_history_set_t), dimension(:), allocatable :: &
         res_history_set
    res_history_set = prc_set%res_history_set
  end function resonant_subprocess_set_get_resonance_history_set

  elemental module function get_libname_res (proc_id) result (libname)
    type(string_t), intent(in) :: proc_id
    type(string_t) :: libname
    libname = proc_id // "_R"
  end function get_libname_res

  module subroutine spawn_resonant_subprocess_libraries &
       (libname, local, global, libname_res)
    type(string_t), intent(in) :: libname
    type(rt_data_t), intent(inout), target :: local
    type(rt_data_t), intent(inout), target :: global
    type(string_t), dimension(:), allocatable, intent(inout) :: libname_res
    type(process_library_t), pointer :: lib
    type(string_t), dimension(:), allocatable :: process_id_res
    type(process_t), pointer :: process
    type(resonance_history_set_t) :: res_history_set
    type(process_component_def_t), pointer :: process_component_def
    logical :: phs_only_saved, exist
    integer :: i_proc, i_component
    lib => global%prclib_stack%get_library_ptr (libname)
    call lib%get_process_id_req_resonant (process_id_res)
    if (size (process_id_res) > 0) then
       call msg_message ("Creating resonant-subprocess libraries &
            &for library '" // char (libname) // "'")
       libname_res = get_libname_res (process_id_res)
       phs_only_saved = local%var_list%get_lval (var_str ("?phs_only"))
       call local%var_list%set_log &
            (var_str ("?phs_only"), .true., is_known=.true.)
       do i_proc = 1, size (process_id_res)
          associate (proc_id => process_id_res (i_proc))
            call msg_message ("Process '" // char (proc_id) // "': &
                 &constructing phase space for resonance structure")
            call integrate_process (proc_id, local, global)
            process => global%process_stack%get_process_ptr (proc_id)
            call create_library (libname_res(i_proc), global, exist)
            if (.not. exist) then
               do i_component = 1, process%get_n_components ()
                  call process%extract_resonance_history_set &
                       (res_history_set, i_component = i_component)
                  process_component_def &
                       => process%get_component_def_ptr (i_component)
                  call add_to_library (libname_res(i_proc), &
                       res_history_set, &
                       process_component_def%get_prt_spec_in (), &
                       process_component_def%get_prt_spec_out (), &
                       global)
               end do
               call msg_message ("Process library '" &
                    // char (libname_res(i_proc)) &
                    // "': created")
            end if
            call global%update_prclib (lib)
          end associate
       end do
       call local%var_list%set_log &
            (var_str ("?phs_only"), phs_only_saved, is_known=.true.)
    end if
  end subroutine spawn_resonant_subprocess_libraries

  module subroutine resonant_subprocess_set_create_library (prc_set, &
       libname, global, exist)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    type(string_t), intent(in) :: libname
    type(rt_data_t), intent(inout), target :: global
    logical, intent(out) :: exist
    prc_set%libname = libname
    call create_library (prc_set%libname, global, exist)
  end subroutine resonant_subprocess_set_create_library

  module subroutine resonant_subprocess_set_add_to_library (prc_set, &
       i_component, prt_in, prt_out, global)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    integer, intent(in) :: i_component
    type(prt_spec_t), dimension(:), intent(in) :: prt_in
    type(prt_spec_t), dimension(:), intent(in) :: prt_out
    type(rt_data_t), intent(inout), target :: global
    call add_to_library (prc_set%libname, &
         prc_set%res_history_set(i_component), &
         prt_in, prt_out, global)
  end subroutine resonant_subprocess_set_add_to_library

  module subroutine resonant_subprocess_set_freeze_library (prc_set, global)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    type(rt_data_t), intent(inout), target :: global
    type(prclib_entry_t), pointer :: lib_entry
    type(process_library_t), pointer :: lib
    lib => global%prclib_stack%get_library_ptr (prc_set%libname)
    call lib%get_process_id_list (prc_set%proc_id)
    prc_set%lib_active = .true.
  end subroutine resonant_subprocess_set_freeze_library

  subroutine create_library (libname, global, exist)
    type(string_t), intent(in) :: libname
    type(rt_data_t), intent(inout), target :: global
    logical, intent(out) :: exist
    type(prclib_entry_t), pointer :: lib_entry
    type(process_library_t), pointer :: lib
    type(resonance_history_t) :: res_history
    type(string_t), dimension(:), allocatable :: proc_id
    type(restricted_process_configuration_t) :: prc_config
    integer :: i
    lib => global%prclib_stack%get_library_ptr (libname)
    exist = associated (lib)
    if (.not. exist) then
       call msg_message ("Creating library for resonant subprocesses '" &
            // char (libname) // "'")
       allocate (lib_entry)
       call lib_entry%init (libname)
       lib => lib_entry%process_library_t
       call global%add_prclib (lib_entry)
    else
       call msg_message ("Using library for resonant subprocesses '" &
            // char (libname) // "'")
       call global%update_prclib (lib)
    end if
  end subroutine create_library

  subroutine add_to_library (libname, res_history_set, prt_in, prt_out, global)
    type(string_t), intent(in) :: libname
    type(resonance_history_set_t), intent(in) :: res_history_set
    type(prt_spec_t), dimension(:), intent(in) :: prt_in
    type(prt_spec_t), dimension(:), intent(in) :: prt_out
    type(rt_data_t), intent(inout), target :: global
    type(prclib_entry_t), pointer :: lib_entry
    type(process_library_t), pointer :: lib
    type(resonance_history_t) :: res_history
    type(string_t), dimension(:), allocatable :: proc_id
    type(restricted_process_configuration_t) :: prc_config
    integer :: n0, i
    lib => global%prclib_stack%get_library_ptr (libname)
    if (associated (lib)) then
       n0 = lib%get_n_processes ()
       allocate (proc_id (res_history_set%get_n_history ()))
       do i = 1, size (proc_id)
          proc_id(i) = libname // str (n0 + i)
          res_history = res_history_set%get_history(i)
          call prc_config%init_resonant_process (proc_id(i), &
               prt_in, prt_out, &
               res_history, &
               global%model, global%var_list)
          call msg_message ("Resonant subprocess #" &
               // char (str(n0+i)) // ": " &
               // char (res_history%as_omega_string (size (prt_in))))
          call prc_config%record (global)
          if (signal_is_pending ())  return
       end do
    else
       call msg_bug ("Adding subprocesses: library '" &
            // char (libname) // "' not found")
    end if
  end subroutine add_to_library

  module subroutine resonant_subprocess_set_compile_library (prc_set, global)
    class(resonant_subprocess_set_t), intent(in) :: prc_set
    type(rt_data_t), intent(inout), target :: global
    type(process_library_t), pointer :: lib
    lib => global%prclib_stack%get_library_ptr (prc_set%libname)
    if (lib%get_status () < STAT_ACTIVE) then
       call compile_library (prc_set%libname, global)
    end if
  end subroutine resonant_subprocess_set_compile_library

  module function resonant_subprocess_set_is_active (prc_set) result (flag)
    class(resonant_subprocess_set_t), intent(in) :: prc_set
    logical :: flag
    flag = prc_set%lib_active
  end function resonant_subprocess_set_is_active

  module function resonant_subprocess_set_get_n_process (prc_set) result (n)
    class(resonant_subprocess_set_t), intent(in) :: prc_set
    integer :: n
    if (prc_set%lib_active) then
       n = size (prc_set%proc_id)
    else
       n = 0
    end if
  end function resonant_subprocess_set_get_n_process

  module function resonant_subprocess_set_get_libname (prc_set) result (libname)
    class(resonant_subprocess_set_t), intent(in) :: prc_set
    type(string_t) :: libname
    if (prc_set%lib_active) then
       libname = prc_set%libname
    else
       libname = ""
    end if
  end function resonant_subprocess_set_get_libname

  module function resonant_subprocess_set_get_proc_id &
       (prc_set, i) result (proc_id)
    class(resonant_subprocess_set_t), intent(in) :: prc_set
    integer, intent(in) :: i
    type(string_t) :: proc_id
    if (allocated (prc_set%proc_id)) then
       proc_id = prc_set%proc_id(i)
    else
       proc_id = ""
    end if
  end function resonant_subprocess_set_get_proc_id

  module subroutine resonant_subprocess_set_prepare_process_objects &
       (prc_set, local, global)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    type(rt_data_t), intent(inout), target :: local
    type(rt_data_t), intent(inout), optional, target :: global
    type(rt_data_t), pointer :: current
    type(process_library_t), pointer :: lib
    type(string_t) :: phs_method_saved, integration_method_saved
    type(string_t) :: proc_id, libname_cur, libname_res
    integer :: i, n
    if (.not. prc_set%is_active ())  return
    if (present (global)) then
       current => global
    else
       current => local
    end if
    libname_cur = current%prclib%get_name ()
    libname_res = prc_set%get_libname ()
    lib => current%prclib_stack%get_library_ptr (libname_res)
    if (associated (lib))  call current%update_prclib (lib)
    phs_method_saved = local%get_sval (var_str ("$phs_method"))
    integration_method_saved = local%get_sval (var_str ("$integration_method"))
    call local%set_string (var_str ("$phs_method"), &
            var_str ("none"), is_known = .true.)
    call local%set_string (var_str ("$integration_method"), &
            var_str ("none"), is_known = .true.)
    n = prc_set%get_n_process ()
    allocate (prc_set%subprocess (n))
    do i = 1, n
       proc_id = prc_set%get_proc_id (i)
       call prepare_process (prc_set%subprocess(i)%p, proc_id)
       if (signal_is_pending ())  return
    end do
    call local%set_string (var_str ("$phs_method"), &
            phs_method_saved, is_known = .true.)
    call local%set_string (var_str ("$integration_method"), &
            integration_method_saved, is_known = .true.)
    lib => current%prclib_stack%get_library_ptr (libname_cur)
    if (associated (lib))  call current%update_prclib (lib)
  contains
    subroutine prepare_process (process, process_id)
      type(process_t), pointer, intent(out) :: process
      type(string_t), intent(in) :: process_id
      call msg_message ("Simulate: initializing resonant subprocess '" &
               // char (process_id) // "'")
      if (present (global)) then
         call integrate_process (process_id, local, global, &
              init_only = .true.)
      else
         call integrate_process (process_id, local, local_stack = .true., &
              init_only = .true.)
      end if
      process => current%process_stack%get_process_ptr (process_id)
      if (.not. associated (process)) then
         call msg_fatal ("Simulate: resonant subprocess '" &
               // char (process_id) // "' could not be initialized: aborting")
      end if
    end subroutine prepare_process
  end subroutine resonant_subprocess_set_prepare_process_objects

  module subroutine resonant_subprocess_set_prepare_process_instances &
       (prc_set, global)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    type(rt_data_t), intent(in), target :: global
    integer :: i, n
    if (.not. prc_set%is_active ())  return
    n = size (prc_set%subprocess)
    allocate (prc_set%instance (n))
    do i = 1, n
       allocate (prc_set%instance(i)%p)
       call prc_set%instance(i)%p%init (prc_set%subprocess(i)%p)
       call prc_set%instance(i)%p%setup_event_data (global%model)
    end do
  end subroutine resonant_subprocess_set_prepare_process_instances

  module subroutine resonant_subprocess_set_connect_transform (prc_set, evt)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    class(evt_t), intent(in), target :: evt
    select type (evt)
    type is (evt_resonance_t)
       prc_set%evt => evt
       call prc_set%evt%set_subprocess_instances (prc_set%instance)
    class default
       call msg_bug ("Resonant subprocess set: event transform has wrong type")
    end select
  end subroutine resonant_subprocess_set_connect_transform

  module subroutine resonant_subprocess_set_on_shell_limit &
       (prc_set, on_shell_limit)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    real(default), intent(in) :: on_shell_limit
    call prc_set%evt%set_on_shell_limit (on_shell_limit)
  end subroutine resonant_subprocess_set_on_shell_limit

  module subroutine resonant_subprocess_set_on_shell_turnoff &
       (prc_set, on_shell_turnoff)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    real(default), intent(in) :: on_shell_turnoff
    call prc_set%evt%set_on_shell_turnoff (on_shell_turnoff)
  end subroutine resonant_subprocess_set_on_shell_turnoff

  module subroutine resonant_subprocess_set_background_factor &
       (prc_set, background_factor)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    real(default), intent(in) :: background_factor
    call prc_set%evt%set_background_factor (background_factor)
  end subroutine resonant_subprocess_set_background_factor

  module subroutine resonant_subprocess_set_dump_instances &
       (prc_set, unit, testflag)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: i, n, u
    u = given_output_unit (unit)
    write (u, "(A)")  "*** Process instances of resonant subprocesses"
    write (u, *)
    n = size (prc_set%subprocess)
    do i = 1, n
       associate (instance => prc_set%instance(i)%p)
         call instance%write (u, testflag)
         write (u, *)
         write (u, *)
       end associate
    end do
  end subroutine resonant_subprocess_set_dump_instances

  module subroutine resonant_subprocess_set_fill_momenta (prc_set)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    integer :: i, n
    call prc_set%evt%fill_momenta ()
  end subroutine resonant_subprocess_set_fill_momenta

  module subroutine resonant_subprocess_set_determine_on_shell_histories &
       (prc_set, i_component, index_array)
    class(resonant_subprocess_set_t), intent(in) :: prc_set
    integer, intent(in) :: i_component
    integer, dimension(:), allocatable, intent(out) :: index_array
    call prc_set%evt%determine_on_shell_histories (index_array)
  end subroutine resonant_subprocess_set_determine_on_shell_histories

  module subroutine resonant_subprocess_set_evaluate_subprocess &
       (prc_set, index_array)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    integer, dimension(:), intent(in) :: index_array
    call prc_set%evt%evaluate_subprocess (index_array)
  end subroutine resonant_subprocess_set_evaluate_subprocess

  module function resonant_subprocess_set_get_master_sqme &
       (prc_set) result (sqme)
    class(resonant_subprocess_set_t), intent(in) :: prc_set
    real(default) :: sqme
    sqme = prc_set%evt%get_master_sqme ()
  end function resonant_subprocess_set_get_master_sqme

  module subroutine resonant_subprocess_set_get_subprocess_sqme (prc_set, sqme)
    class(resonant_subprocess_set_t), intent(in) :: prc_set
    real(default), dimension(:), intent(inout) :: sqme
    integer :: i
    call prc_set%evt%get_subprocess_sqme (sqme)
  end subroutine resonant_subprocess_set_get_subprocess_sqme

  module subroutine resonant_subprocess_set_compute_probabilities &
       (prc_set, prob_array)
    class(resonant_subprocess_set_t), intent(inout) :: prc_set
    real(default), dimension(:), allocatable, intent(out) :: prob_array
    integer, dimension(:), allocatable :: index_array
    real(default) :: sqme, sqme_sum, sqme_bg
    real(default), dimension(:), allocatable :: sqme_res
    integer :: n
    n = size (prc_set%subprocess)
    allocate (prob_array (0:n), source = 0._default)
    call prc_set%evt%compute_probabilities ()
    call prc_set%evt%get_selector_weights (prob_array)
  end subroutine resonant_subprocess_set_compute_probabilities


end submodule restricted_subprocesses_s

