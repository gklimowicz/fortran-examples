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

submodule (process_configurations) process_configurations_s

  use debug_master, only: debug_on
  use diagnostics
  use io_units
  use physics_defs, only: BORN, NLO_VIRTUAL, NLO_REAL, NLO_DGLAP, &
       NLO_SUBTRACTION, NLO_MISMATCH
  use prc_core_def
  use dispatch_me_methods, only: dispatch_core_def
  use prc_external, only: prc_external_def_t

  implicit none

contains

  module subroutine process_configuration_write (config, unit)
    class(process_configuration_t), intent(in) :: config
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(A)")  "Process configuration:"
    if (associated (config%entry)) then
       call config%entry%write (u)
    else
       write (u, "(1x,3A)")  "ID    = '", char (config%id), "'"
       write (u, "(1x,A,1x,I0)")  "num ID =", config%num_id
       write (u, "(2x,A)")  "[no entry]"
    end if
  end subroutine process_configuration_write

  module subroutine process_configuration_init &
       (config, prc_name, n_in, n_components, model, var_list, &
       nlo_process, negative_sf)
    class(process_configuration_t), intent(out) :: config
    type(string_t), intent(in) :: prc_name
    integer, intent(in) :: n_in
    integer, intent(in) :: n_components
    type(model_t), intent(in), pointer :: model
    type(var_list_t), intent(in) :: var_list
    logical, intent(in), optional :: nlo_process, negative_sf
    logical :: nlo_proc, neg_sf
    logical :: requires_resonances
    if (debug_on) call msg_debug (D_CORE, "process_configuration_init")
    config%id = prc_name
    if (present (nlo_process)) then
       nlo_proc = nlo_process
    else
       nlo_proc = .false.
    end if
    if (present (negative_sf)) then
       neg_sf = negative_sf
    else
       neg_sf = .false.
    end if
    requires_resonances = var_list%get_lval (var_str ("?resonance_history"))

    if (debug_on) call msg_debug (D_CORE, "nlo_process", nlo_proc)
    allocate (config%entry)
    if (var_list%is_known (var_str ("process_num_id"))) then
       config%num_id = &
            var_list%get_ival (var_str ("process_num_id"))
       call config%entry%init (prc_name, &
            model = model, n_in = n_in, n_components = n_components, &
            num_id = config%num_id, &
            nlo_process = nlo_proc, &
            negative_sf = neg_sf, &
            requires_resonances = requires_resonances)
    else
       call config%entry%init (prc_name, &
            model = model, n_in = n_in, n_components = n_components, &
            nlo_process = nlo_proc, &
            negative_sf = neg_sf, &
            requires_resonances = requires_resonances)
    end if
  end subroutine process_configuration_init

  module subroutine process_configuration_setup_component &
       (config, i_component, prt_in, prt_out, model, var_list, &
        nlo_type, can_be_integrated)
    class(process_configuration_t), intent(inout) :: config
    integer, intent(in) :: i_component
    type(prt_spec_t), dimension(:), intent(in) :: prt_in
    type(prt_spec_t), dimension(:), intent(in) :: prt_out
    type(model_t), pointer, intent(in) :: model
    type(var_list_t), intent(in) :: var_list
    integer, intent(in), optional :: nlo_type
    logical, intent(in), optional :: can_be_integrated
    type(string_t), dimension(:), allocatable :: prt_str_in
    type(string_t), dimension(:), allocatable :: prt_str_out
    class(prc_core_def_t), allocatable :: core_def
    type(string_t) :: method
    type(string_t) :: born_me_method
    type(string_t) :: real_tree_me_method
    type(string_t) :: loop_me_method
    type(string_t) :: correlation_me_method
    type(string_t) :: dglap_me_method
    integer :: i
    if (debug_on)  call msg_debug2 &
         (D_CORE, "process_configuration_setup_component")
    allocate (prt_str_in  (size (prt_in)))
    allocate (prt_str_out (size (prt_out)))
    forall (i = 1:size (prt_in))  prt_str_in(i)  = prt_in(i)% get_name ()
    forall (i = 1:size (prt_out)) prt_str_out(i) = prt_out(i)%get_name ()

    method = var_list%get_sval (var_str ("$method"))
    if (present (nlo_type)) then
       select case (nlo_type)
       case (BORN)
          born_me_method = var_list%get_sval (var_str ("$born_me_method"))
          if (born_me_method /= var_str ("")) then
             method = born_me_method
          end if
       case (NLO_VIRTUAL)
          loop_me_method = var_list%get_sval (var_str ("$loop_me_method"))
          if (loop_me_method /= var_str ("")) then
             method = loop_me_method
          end if
       case (NLO_REAL)
          real_tree_me_method = &
               var_list%get_sval (var_str ("$real_tree_me_method"))
          if (real_tree_me_method /= var_str ("")) then
             method = real_tree_me_method
          end if
       case (NLO_DGLAP)
          dglap_me_method = &
               var_list%get_sval (var_str ("$dglap_me_method"))
          if (dglap_me_method /= var_str ("")) then
             method = dglap_me_method
          end if
       case (NLO_SUBTRACTION,NLO_MISMATCH)
          correlation_me_method = &
               var_list%get_sval (var_str ("$correlation_me_method"))
          if (correlation_me_method /= var_str ("")) then
             method = correlation_me_method
          end if
       case default
       end select
    end if
    call dispatch_core_def (core_def, prt_str_in, prt_str_out, &
         model, var_list, config%id, nlo_type, method)
    select type (core_def)
    class is (prc_external_def_t)
       if (present (can_be_integrated)) then
          call core_def%set_active_writer (can_be_integrated)
       else
          call msg_fatal ("Cannot decide if external core is integrated!")
       end if
    end select

    if (debug_on)  call msg_debug2 &
         (D_CORE, "import_component with method ", method)
    call config%entry%import_component (i_component, &
         n_out = size (prt_out), &
         prt_in = prt_in, &
         prt_out = prt_out, &
         method = method, &
         variant = core_def, &
         nlo_type = nlo_type, &
         can_be_integrated = can_be_integrated)
  end subroutine process_configuration_setup_component

  module subroutine process_configuration_set_fixed_emitter (config, i, emitter)
     class(process_configuration_t), intent(inout) :: config
     integer, intent(in) :: i, emitter
     call config%entry%set_fixed_emitter (i, emitter)
  end subroutine process_configuration_set_fixed_emitter

  module subroutine process_configuration_set_coupling_powers &
       (config, alpha_power, alphas_power)
    class(process_configuration_t), intent(inout) :: config
    integer, intent(in) :: alpha_power, alphas_power
    call config%entry%set_coupling_powers (alpha_power, alphas_power)
  end subroutine process_configuration_set_coupling_powers

  module subroutine process_configuration_set_component_associations &
         (config, i_list, remnant, use_real_finite, mismatch)
    class(process_configuration_t), intent(inout) :: config
    integer, dimension(:), intent(in) :: i_list
    logical, intent(in) :: remnant, use_real_finite, mismatch
    integer :: i_component
    do i_component = 1, config%entry%get_n_components ()
       if (any (i_list == i_component)) then
          call config%entry%set_associated_components (i_component, &
               i_list, remnant, use_real_finite, mismatch)
       end if
    end do
  end subroutine process_configuration_set_component_associations

  module subroutine process_configuration_record (config, global)
    class(process_configuration_t), intent(inout) :: config
    type(rt_data_t), intent(inout) :: global
    if (associated (global%prclib)) then
       call global%prclib%open ()
       call global%prclib%append (config%entry)
       if (config%num_id /= 0) then
          write (msg_buffer, "(5A,I0,A)") "Process library '", &
               char (global%prclib%get_name ()), &
               "': recorded process '", char (config%id), "' (", &
               config%num_id, ")"
       else
          write (msg_buffer, "(5A)") "Process library '", &
               char (global%prclib%get_name ()), &
               "': recorded process '", char (config%id), "'"
       end if
       call msg_message ()
    else
       call msg_fatal ("Recording process '" // char (config%id) &
            // "': active process library undefined")
    end if
  end subroutine process_configuration_record


end submodule process_configurations_s

