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

submodule (dispatch_mci) dispatch_mci_s

  use diagnostics
  use os_interface
  use mci_none
  use mci_midpoint
  use mci_vamp
  use mci_vamp2

  implicit none

  character(*), parameter :: ALLOWED_IN_DIRNAME = &
       "abcdefghijklmnopqrstuvwxyz&
       &ABCDEFGHIJKLMNOPQRSTUVWXYZ&
       &1234567890&
       &.,_-+="

contains

  module subroutine dispatch_mci_setup (mci, var_list, process_id, is_nlo)
    class(mci_t), allocatable, intent(out) :: mci
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_id
    logical, intent(in), optional :: is_nlo
    type(string_t) :: run_id
    type(string_t) :: integration_method
    type(grid_parameters_t) :: grid_par
    type(history_parameters_t) :: history_par
    type(mci_vamp2_config_t) :: mci_vamp2_config
    integer :: grid_checkpoint
    logical :: rebuild_grids, check_grid_file, negative_weights, verbose
    logical :: dispatch_nlo, binary_grid_format
    type(string_t) :: grid_path, parallel_method
    dispatch_nlo = .false.; if (present (is_nlo)) dispatch_nlo = is_nlo
    integration_method = &
         var_list%get_sval (var_str ("$integration_method"))
    select case (char (integration_method))
    case ("none")
       allocate (mci_none_t :: mci)
    case ("midpoint")
       allocate (mci_midpoint_t :: mci)
    case ("vamp", "default")
       call unpack_options_vamp ()
       allocate (mci_vamp_t :: mci)
       select type (mci)
       type is (mci_vamp_t)
          call mci%set_grid_parameters (grid_par)
          if (run_id /= "") then
             call mci%set_grid_filename (process_id, run_id)
          else
             call mci%set_grid_filename (process_id)
          end if
          grid_path = var_list%get_sval (var_str ("$integrate_workspace"))
          if (grid_path /= "") then
             call setup_grid_path (grid_path)
             call mci%prepend_grid_path (grid_path)
          end if
          call mci%set_history_parameters (history_par)
          call mci%set_rebuild_flag (rebuild_grids, check_grid_file)
          mci%negative_weights = negative_weights
          mci%verbose = verbose
       end select
    case ("vamp2")
       call unpack_options_vamp2 ()
       allocate (mci_vamp2_t :: mci)
       select type (mci)
       type is (mci_vamp2_t)
          call mci%set_config (mci_vamp2_config)
          if (run_id /= "") then
             call mci%set_grid_filename (process_id, run_id)
          else
             call mci%set_grid_filename (process_id)
          end if
          grid_path = var_list%get_sval (var_str ("$integrate_workspace"))
          if (grid_path /= "") then
             call setup_grid_path (grid_path)
             call mci%prepend_grid_path (grid_path)
          end if
          call mci%set_rebuild_flag (rebuild_grids, check_grid_file)
          mci%negative_weights = negative_weights
          mci%verbose = verbose
          mci%grid_checkpoint = grid_checkpoint
          mci%binary_grid_format = binary_grid_format
          mci%parallel_method = parallel_method
       end select
    case default
       call msg_fatal ("Integrator '" &
            // char (integration_method) // "' not implemented")
    end select
  contains
      subroutine unpack_options_vamp ()
        grid_par%threshold_calls = &
             var_list%get_ival (var_str ("threshold_calls"))
        grid_par%min_calls_per_channel = &
             var_list%get_ival (var_str ("min_calls_per_channel"))
        grid_par%min_calls_per_bin = &
             var_list%get_ival (var_str ("min_calls_per_bin"))
        grid_par%min_bins = &
             var_list%get_ival (var_str ("min_bins"))
        grid_par%max_bins = &
             var_list%get_ival (var_str ("max_bins"))
        grid_par%stratified = &
             var_list%get_lval (var_str ("?stratified"))
        select case (char (var_list%get_sval (var_str ("$phs_method"))))
        case ("rambo")
           grid_par%use_vamp_equivalences = .false.
        case default
           grid_par%use_vamp_equivalences = &
           var_list%get_lval (var_str ("?use_vamp_equivalences"))
        end select
        grid_par%channel_weights_power = &
             var_list%get_rval (var_str ("channel_weights_power"))
        grid_par%accuracy_goal = &
             var_list%get_rval (var_str ("accuracy_goal"))
        grid_par%error_goal = &
             var_list%get_rval (var_str ("error_goal"))
        grid_par%rel_error_goal = &
             var_list%get_rval (var_str ("relative_error_goal"))
        history_par%global = &
             var_list%get_lval (var_str ("?vamp_history_global"))
        history_par%global_verbose = &
             var_list%get_lval (var_str ("?vamp_history_global_verbose"))
        history_par%channel = &
             var_list%get_lval (var_str ("?vamp_history_channels"))
        history_par%channel_verbose = &
             var_list%get_lval (var_str ("?vamp_history_channels_verbose"))
        verbose = &
             var_list%get_lval (var_str ("?vamp_verbose"))
        check_grid_file = &
             var_list%get_lval (var_str ("?check_grid_file"))
        run_id = &
             var_list%get_sval (var_str ("$run_id"))
        rebuild_grids = &
             var_list%get_lval (var_str ("?rebuild_grids"))
        negative_weights = &
             var_list%get_lval (var_str ("?negative_weights")) .or. dispatch_nlo
      end subroutine unpack_options_vamp

      subroutine unpack_options_vamp2 ()
        mci_vamp2_config%n_bins_max = &
             var_list%get_ival (var_str ("max_bins"))
        mci_vamp2_config%n_calls_min_per_channel = &
             var_list%get_ival (var_str ("min_calls_per_channel"))
        mci_vamp2_config%n_calls_threshold = &
             var_list%get_ival (var_str ("threshold_calls"))
        mci_vamp2_config%beta = &
             var_list%get_rval (var_str ("channel_weights_power"))
        mci_vamp2_config%stratified = &
             var_list%get_lval (var_str ("?stratified"))
        select case (char (var_list%get_sval (var_str ("$phs_method"))))
        case ("rambo")
           mci_vamp2_config%equivalences = .false.
        case default
           mci_vamp2_config%equivalences = &
              var_list%get_lval (var_str ("?use_vamp_equivalences"))
        end select
        mci_vamp2_config%accuracy_goal = &
             var_list%get_rval (var_str ("accuracy_goal"))
        mci_vamp2_config%error_goal = &
             var_list%get_rval (var_str ("error_goal"))
        mci_vamp2_config%rel_error_goal = &
             var_list%get_rval (var_str ("relative_error_goal"))
        verbose = &
             var_list%get_lval (var_str ("?vamp_verbose"))
        check_grid_file = &
             var_list%get_lval (var_str ("?check_grid_file"))
        run_id = &
             var_list%get_sval (var_str ("$run_id"))
        rebuild_grids = &
             var_list%get_lval (var_str ("?rebuild_grids"))
        negative_weights = &
             var_list%get_lval (var_str ("?negative_weights")) .or. dispatch_nlo
        grid_checkpoint = &
             var_list%get_ival (var_str ("vamp_grid_checkpoint"))
        select case (char (var_list%get_sval (var_str ("$vamp_grid_format"))))
        case ("binary","Binary","BINARY")
           binary_grid_format = .true.
        case ("ascii","Ascii","ASCII")
           binary_grid_format = .false.
        case default
           binary_grid_format = .false.
        end select
        select case (char (var_list%get_sval (var_str ("$vamp_parallel_method"))))
        case ("simple","Simple","SIMPLE")
           parallel_method = var_str ("simple")
        case ("load","Load","LOAD")
           parallel_method = var_str ("load")
        case default
           parallel_method = var_str ("simple")
        end select
      end subroutine unpack_options_vamp2

  end subroutine dispatch_mci_setup

  module subroutine setup_grid_path (grid_path)
    type(string_t), intent(in) :: grid_path
    if (verify (grid_path, ALLOWED_IN_DIRNAME) == 0) then
       call msg_message ("Integrator: preparing VAMP grid directory '" &
            // char (grid_path) // "'")
       call os_system_call ("mkdir -p '" // grid_path // "'")
    else
       call msg_fatal ("Integrator: VAMP grid_path '" &
            // char (grid_path) // "' contains illegal characters")
    end if
  end subroutine setup_grid_path


end submodule dispatch_mci_s

