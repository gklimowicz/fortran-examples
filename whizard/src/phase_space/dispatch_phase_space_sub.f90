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

submodule (dispatch_phase_space) dispatch_phase_space_s

  use io_units, only: free_unit
  use diagnostics
  use phs_none
  use phs_single
  use phs_rambo
  use phs_wood
  use phs_fks

  implicit none

contains

  module subroutine dispatch_phs (phs, var_list, os_data, process_id, &
         mapping_defaults, phs_par, phs_method_in)
    class(phs_config_t), allocatable, intent(inout) :: phs
    type(var_list_t), intent(in) :: var_list
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in) :: process_id
    type(mapping_defaults_t), intent(in), optional :: mapping_defaults
    type(phs_parameters_t), intent(in), optional :: phs_par
    type(string_t), intent(in), optional :: phs_method_in
    type(string_t) :: phs_method, phs_file, run_id
    logical :: use_equivalences, vis_channels, fatal_beam_decay
    integer :: u_phs
    logical :: exist
    if (present (phs_method_in)) then
       phs_method = phs_method_in
    else
       phs_method = &
            var_list%get_sval (var_str ("$phs_method"))
    end if
    phs_file = &
         var_list%get_sval (var_str ("$phs_file"))
    use_equivalences = &
         var_list%get_lval (var_str ("?use_vamp_equivalences"))
    vis_channels = &
         var_list%get_lval (var_str ("?vis_channels"))
    fatal_beam_decay = &
         var_list%get_lval (var_str ("?fatal_beam_decay"))
    run_id = &
         var_list%get_sval (var_str ("$run_id"))
    select case (char (phs_method))
    case ("none")
       allocate (phs_none_config_t :: phs)
    case ("single")
       allocate (phs_single_config_t :: phs)
       if (vis_channels) then
          call msg_warning ("Visualizing phase space channels not " // &
               "available for method 'single'.")
       end if
    case ("rambo")
       allocate (phs_rambo_config_t :: phs)
       if (vis_channels) &
          call msg_warning ("Visualizing phase space channels not " // &
              "available for method 'rambo'.")
    case ("fks")
       allocate (phs_fks_config_t :: phs)
       if (use_equivalences) then
          select type (phs)
          type is (phs_fks_config_t)
             call phs%enable_equivalences ()
          end select
       end if
    case ("wood", "default", "fast_wood")
       call dispatch_wood ()
    case default
       call msg_fatal ("Phase space: parameterization method '" &
            // char (phs_method) // "' not implemented")
    end select
  contains
    subroutine dispatch_wood ()
      allocate (phs_wood_config_t :: phs)
      select type (phs)
      type is (phs_wood_config_t)
         if (phs_file /= "") then
            inquire (file = char (phs_file), exist = exist)
            if (exist) then
               call msg_message ("Phase space: reading configuration from '" &
                    // char (phs_file) // "'")
               u_phs = free_unit ()
               open (u_phs, file = char (phs_file), &
                    action = "read", status = "old")
               call phs%set_input (u_phs)
            else
               call msg_fatal ("Phase space: configuration file '" &
                    // char (phs_file) // "' not found")
            end if
         end if
         if (present (phs_par)) &
              call phs%set_parameters (phs_par)
         if (use_equivalences) &
              call phs%enable_equivalences ()
         if (present (mapping_defaults)) &
              call phs%set_mapping_defaults (mapping_defaults)
         if (phs_method == "fast_wood") phs%use_cascades2 = .true.
         phs%vis_channels = vis_channels
         phs%fatal_beam_decay = fatal_beam_decay
         phs%os_data = os_data
         phs%run_id = run_id
      end select
    end subroutine dispatch_wood

  end subroutine dispatch_phs

  module subroutine dispatch_sf_channels (sf_channel, sf_string, sf_prop, &
       coll, var_list, sqrts, beam_structure)
    type(sf_channel_t), dimension(:), allocatable, intent(out) :: sf_channel
    type(string_t), intent(out) :: sf_string
    type(sf_prop_t), intent(in) :: sf_prop
    type(phs_channel_collection_t), intent(in) :: coll
    type(var_list_t), intent(in) :: var_list
    real(default), intent(in) :: sqrts
    type(beam_structure_t), intent(in) :: beam_structure
    type(beam_structure_t) :: beam_structure_tmp
    class(channel_prop_t), allocatable :: prop
    integer :: n_strfun, n_sf_channel, i
    logical :: sf_allow_s_mapping, circe1_map, circe1_generate
    logical :: s_mapping_enable, endpoint_mapping, power_mapping
    logical :: single_parameter
    integer, dimension(:), allocatable :: s_mapping, single_mapping
    real(default) :: s_mapping_power
    real(default) :: circe1_mapping_slope, endpoint_mapping_slope
    real(default) :: power_mapping_eps
    beam_structure_tmp = beam_structure
    call beam_structure_tmp%expand (strfun_mode)
    n_strfun = beam_structure_tmp%get_n_record ()
    sf_string = beam_structure_tmp%to_string (sf_only = .true.)
    sf_allow_s_mapping = &
         var_list%get_lval (var_str ("?sf_allow_s_mapping"))
    circe1_generate = &
         var_list%get_lval (var_str ("?circe1_generate"))
    circe1_map = &
         var_list%get_lval (var_str ("?circe1_map"))
    circe1_mapping_slope = &
         var_list%get_rval (var_str ("circe1_mapping_slope"))
    s_mapping_enable = .false.
    s_mapping_power = 1
    endpoint_mapping = .false.
    endpoint_mapping_slope = 1
    power_mapping = .false.
    single_parameter = .false.
    select case (char (sf_string))
    case ("", "[any particles]")
    case ("pdf_builtin, none", &
         "pdf_builtin_photon, none", &
         "none, pdf_builtin", &
         "none, pdf_builtin_photon", &
         "lhapdf, none", &
         "lhapdf_photon, none", &
         "none, lhapdf", &
         "none, lhapdf_photon")
         single_parameter = .true.
    case ("pdf_builtin, none => none, pdf_builtin", &
          "pdf_builtin, none => none, pdf_builtin_photon", &
          "pdf_builtin_photon, none => none, pdf_builtin", &
          "pdf_builtin_photon, none => none, pdf_builtin_photon", &
          "lhapdf, none => none, lhapdf", &
          "lhapdf, none => none, lhapdf_photon", &
          "lhapdf_photon, none => none, lhapdf", &
          "lhapdf_photon, none => none, lhapdf_photon")
       allocate (s_mapping (2), source = [1, 2])
       s_mapping_enable = .true.
       s_mapping_power = 2
    case ("pdf_builtin, none => none, pdf_builtin => epa, none => none, epa", &
          "pdf_builtin, none => none, pdf_builtin => ewa, none => none, ewa", &
          "pdf_builtin, none => none, pdf_builtin => ewa, none => none, epa", &
          "pdf_builtin, none => none, pdf_builtin => epa, none => none, ewa")
       allocate (s_mapping (2), source = [1, 2])
       s_mapping_enable = .true.
       s_mapping_power = 2
    case ("isr, none", &
         "none, isr")
       allocate (single_mapping (1), source = [1])
       single_parameter = .true.
    case ("isr, none => none, isr")
       allocate (s_mapping (2), source = [1, 2])
       power_mapping = .true.
       power_mapping_eps = minval (sf_prop%isr_eps)
    case ("isr, none => none, isr => epa, none => none, epa", &
          "isr, none => none, isr => ewa, none => none, ewa", &
          "isr, none => none, isr => ewa, none => none, epa", &
          "isr, none => none, isr => epa, none => none, ewa")
       allocate (s_mapping (2), source = [1, 2])
       power_mapping = .true.
       power_mapping_eps = minval (sf_prop%isr_eps)
    case ("circe1 => isr, none => none, isr => epa, none => none, epa", &
          "circe1 => isr, none => none, isr => ewa, none => none, ewa", &
          "circe1 => isr, none => none, isr => ewa, none => none, epa", &
          "circe1 => isr, none => none, isr => epa, none => none, ewa")
       if (circe1_generate) then
          allocate (s_mapping (2), source = [2, 3])
       else
          allocate (s_mapping (3), source = [1, 2, 3])
          endpoint_mapping = .true.
          endpoint_mapping_slope = circe1_mapping_slope
       end if
       power_mapping = .true.
       power_mapping_eps = minval (sf_prop%isr_eps)
    case ("pdf_builtin, none => none, isr", &
         "pdf_builtin_photon, none => none, isr", &
         "lhapdf, none => none, isr", &
         "lhapdf_photon, none => none, isr")
       allocate (single_mapping (1), source = [2])
    case ("isr, none => none, pdf_builtin", &
         "isr, none => none, pdf_builtin_photon", &
         "isr, none => none, lhapdf", &
         "isr, none => none, lhapdf_photon")
       allocate (single_mapping (1), source = [1])
    case ("epa, none", &
          "none, epa")
       allocate (single_mapping (1), source = [1])
       single_parameter = .true.
    case ("epa, none => none, epa")
       allocate (single_mapping (2), source = [1, 2])
    case ("epa, none => none, isr", &
         "isr, none => none, epa", &
         "ewa, none => none, isr", &
         "isr, none => none, ewa")
       allocate (single_mapping (2), source = [1, 2])
    case ("pdf_builtin, none => none, epa", &
         "pdf_builtin_photon, none => none, epa", &
         "lhapdf, none => none, epa", &
         "lhapdf_photon, none => none, epa")
       allocate (single_mapping (1), source = [2])
    case ("pdf_builtin, none => none, ewa", &
         "pdf_builtin_photon, none => none, ewa", &
         "lhapdf, none => none, ewa", &
         "lhapdf_photon, none => none, ewa")
       allocate (single_mapping (1), source = [2])
    case ("epa, none => none, pdf_builtin", &
         "epa, none => none, pdf_builtin_photon", &
         "epa, none => none, lhapdf", &
         "epa, none => none, lhapdf_photon")
       allocate (single_mapping (1), source = [1])
    case ("ewa, none => none, pdf_builtin", &
         "ewa, none => none, pdf_builtin_photon", &
         "ewa, none => none, lhapdf", &
         "ewa, none => none, lhapdf_photon")
       allocate (single_mapping (1), source = [1])
    case ("ewa, none", &
          "none, ewa")
       allocate (single_mapping (1), source = [1])
       single_parameter = .true.
    case ("ewa, none => none, ewa")
       allocate (single_mapping (2), source = [1, 2])
    case ("energy_scan, none => none, energy_scan")
       allocate (s_mapping (2), source = [1, 2])
    case ("sf_test_1, none => none, sf_test_1")
       allocate (s_mapping (2), source = [1, 2])
    case ("circe1")
       if (circe1_generate) then
          !!! no mapping
       else if (circe1_map) then
          allocate (s_mapping (1), source = [1])
          endpoint_mapping = .true.
          endpoint_mapping_slope = circe1_mapping_slope
       else
          allocate (s_mapping (1), source = [1])
          s_mapping_enable = .true.
       end if
    case ("circe1 => isr, none => none, isr")
       if (circe1_generate) then
          allocate (s_mapping (2), source = [2, 3])
       else
          allocate (s_mapping (3), source = [1, 2, 3])
          endpoint_mapping = .true.
          endpoint_mapping_slope = circe1_mapping_slope
       end if
       power_mapping = .true.
       power_mapping_eps = minval (sf_prop%isr_eps)
    case ("circe1 => isr, none", &
         "circe1 => none, isr")
       allocate (single_mapping (1), source = [2])
    case ("circe1 => epa, none => none, epa")
       if (circe1_generate) then
          allocate (single_mapping (2), source = [2, 3])
       else
          call msg_fatal ("CIRCE/EPA: supported with ?circe1_generate=true &
               &only")
       end if
    case ("circe1 => ewa, none => none, ewa")
       if (circe1_generate) then
          allocate (single_mapping (2), source = [2, 3])
       else
          call msg_fatal ("CIRCE/EWA: supported with ?circe1_generate=true &
               &only")
       end if
    case ("circe1 => epa, none", &
         "circe1 => none, epa")
       if (circe1_generate) then
          allocate (single_mapping (1), source = [2])
       else
          call msg_fatal ("CIRCE/EPA: supported with ?circe1_generate=true &
               &only")
       end if
    case ("circe1 => epa, none => none, isr", &
         "circe1 => isr, none => none, epa", &
         "circe1 => ewa, none => none, isr", &
         "circe1 => isr, none => none, ewa")
       if (circe1_generate) then
          allocate (single_mapping (2), source = [2, 3])
       else
          call msg_fatal ("CIRCE/EPA: supported with ?circe1_generate=true &
               &only")
       end if
    case ("circe2", &
         "gaussian", &
         "beam_events")
       !!! no mapping
    case ("circe2 => isr, none => none, isr", &
       "gaussian => isr, none => none, isr", &
       "beam_events => isr, none => none, isr")
       allocate (s_mapping (2), source = [2, 3])
       power_mapping = .true.
       power_mapping_eps = minval (sf_prop%isr_eps)
    case ("circe2 => isr, none", &
         "circe2 => none, isr", &
         "gaussian => isr, none", &
         "gaussian => none, isr", &
         "beam_events => isr, none", &
         "beam_events => none, isr")
       allocate (single_mapping (1), source = [2])
    case ("circe2 => epa, none => none, epa", &
         "gaussian => epa, none => none, epa", &
         "beam_events => epa, none => none, epa")
       allocate (single_mapping (2), source = [2, 3])
    case ("circe2 => epa, none", &
         "circe2 => none, epa", &
         "circe2 => ewa, none", &
         "circe2 => none, ewa", &
         "gaussian => epa, none", &
         "gaussian => none, epa", &
         "gaussian => ewa, none", &
         "gaussian => none, ewa", &
         "beam_events => epa, none", &
         "beam_events => none, epa", &
         "beam_events => ewa, none", &
         "beam_events => none, ewa")
       allocate (single_mapping (1), source = [2])
    case ("circe2 => epa, none => none, isr", &
         "circe2 => isr, none => none, epa", &
         "circe2 => ewa, none => none, isr", &
         "circe2 => isr, none => none, ewa", &
         "gaussian => epa, none => none, isr", &
         "gaussian => isr, none => none, epa", &
         "gaussian => ewa, none => none, isr", &
         "gaussian => isr, none => none, ewa", &
         "beam_events => epa, none => none, isr", &
         "beam_events => isr, none => none, epa", &
         "beam_events => ewa, none => none, isr", &
         "beam_events => isr, none => none, ewa")
       allocate (single_mapping (2), source = [2, 3])
    case ("energy_scan")
    case default
       call msg_fatal ("Beam structure: " &
            // char (sf_string) // " not supported")
    end select
    if (sf_allow_s_mapping .and. coll%n > 0) then
       n_sf_channel = coll%n
       allocate (sf_channel (n_sf_channel))
       do i = 1, n_sf_channel
          call sf_channel(i)%init (n_strfun)
          if (allocated (single_mapping)) then
             call sf_channel(i)%activate_mapping (single_mapping)
          end if
          if (allocated (prop))  deallocate (prop)
          call coll%get_entry (i, prop)
          if (allocated (prop)) then
             if (endpoint_mapping .and. power_mapping) then
                select type (prop)
                type is (resonance_t)
                   call sf_channel(i)%set_eir_mapping (s_mapping, &
                        a = endpoint_mapping_slope, eps = power_mapping_eps, &
                        m = prop%mass / sqrts, w = prop%width / sqrts)
                type is (on_shell_t)
                   call sf_channel(i)%set_eio_mapping (s_mapping, &
                        a = endpoint_mapping_slope, eps = power_mapping_eps, &
                        m = prop%mass / sqrts)
                end select
             else if (endpoint_mapping) then
                select type (prop)
                type is (resonance_t)
                   call sf_channel(i)%set_epr_mapping (s_mapping, &
                        a = endpoint_mapping_slope, &
                        m = prop%mass / sqrts, w = prop%width / sqrts)
                type is (on_shell_t)
                   call sf_channel(i)%set_epo_mapping (s_mapping, &
                        a = endpoint_mapping_slope, &
                        m = prop%mass / sqrts)
                end select
             else if (power_mapping) then
                select type (prop)
                type is (resonance_t)
                   call sf_channel(i)%set_ipr_mapping (s_mapping, &
                        eps = power_mapping_eps, &
                        m = prop%mass / sqrts, w = prop%width / sqrts)
                type is (on_shell_t)
                   call sf_channel(i)%set_ipo_mapping (s_mapping, &
                        eps = power_mapping_eps, &
                        m = prop%mass / sqrts)
                end select
             else if (allocated (s_mapping)) then
                select type (prop)
                type is (resonance_t)
                   call sf_channel(i)%set_res_mapping (s_mapping, &
                        m = prop%mass / sqrts, w = prop%width / sqrts, &
                        single = single_parameter)
                type is (on_shell_t)
                   call sf_channel(i)%set_os_mapping (s_mapping, &
                        m = prop%mass / sqrts, &
                        single = single_parameter)
                end select
             else if (allocated (single_mapping)) then
                select type (prop)
                type is (resonance_t)
                   call sf_channel(i)%set_res_mapping (single_mapping, &
                           m = prop%mass / sqrts, w = prop%width / sqrts, &
                           single = single_parameter)
                type is (on_shell_t)
                   call sf_channel(i)%set_os_mapping (single_mapping, &
                        m = prop%mass / sqrts, &
                        single = single_parameter)
                end select
             end if
          else if (endpoint_mapping .and. power_mapping) then
             call sf_channel(i)%set_ei_mapping (s_mapping, &
                  a = endpoint_mapping_slope, eps = power_mapping_eps)
          else if (endpoint_mapping .and. .not. allocated (single_mapping)) then
             call sf_channel(i)%set_ep_mapping (s_mapping, &
                  a = endpoint_mapping_slope)
          else if (power_mapping .and. .not. allocated (single_mapping)) then
             call sf_channel(i)%set_ip_mapping (s_mapping, &
                  eps = power_mapping_eps)
          else if (s_mapping_enable .and. .not. allocated (single_mapping)) then
             call sf_channel(i)%set_s_mapping (s_mapping, &
                  power = s_mapping_power)
          end if
       end do
    else if (sf_allow_s_mapping) then
       allocate (sf_channel (1))
       call sf_channel(1)%init (n_strfun)
       if (allocated (single_mapping)) then
          call sf_channel(1)%activate_mapping (single_mapping)
       else if (endpoint_mapping .and. power_mapping) then
          call sf_channel(i)%set_ei_mapping (s_mapping, &
               a = endpoint_mapping_slope, eps = power_mapping_eps)
       else if (endpoint_mapping) then
          call sf_channel(1)%set_ep_mapping (s_mapping, &
                  a = endpoint_mapping_slope)
       else if (power_mapping) then
          call sf_channel(1)%set_ip_mapping (s_mapping, &
                  eps = power_mapping_eps)
       else if (s_mapping_enable) then
          call sf_channel(1)%set_s_mapping (s_mapping, &
               power = s_mapping_power)
       end if
    else
       allocate (sf_channel (1))
       call sf_channel(1)%init (n_strfun)
       if (allocated (single_mapping)) then
          call sf_channel(1)%activate_mapping (single_mapping)
       end if
    end if
  end subroutine dispatch_sf_channels


end submodule dispatch_phase_space_s

