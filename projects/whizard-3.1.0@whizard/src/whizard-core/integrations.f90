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

module integrations

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use diagnostics
  use prc_core
  use process
  use instances
  use process_stacks
  use iterations
  use rt_data
  use nlo_data

  implicit none
  private

  public :: integration_t
  public :: integrate_process

  type :: integration_t
    private
    type(string_t) :: process_id
    type(string_t) :: run_id
    type(process_t), pointer :: process => null ()
    logical :: rebuild_phs = .false.
    logical :: ignore_phs_mismatch = .false.
    logical :: phs_only = .false.
    logical :: process_has_me = .true.
    integer :: n_calls_test = 0
    logical :: vis_history = .true.
    type(string_t) :: history_filename
    type(string_t) :: log_filename
    type(helicity_selection_t) :: helicity_selection
    logical :: use_color_factors = .false.
    logical :: has_beam_pol = .false.
    logical :: combined_integration = .false.
    type(iteration_multipliers_t) :: iteration_multipliers
    type(nlo_settings_t) :: nlo_settings
   contains
     procedure :: create_process => integration_create_process
     procedure :: init_process => integration_init_process
     procedure :: setup_process => integration_setup_process
     procedure :: evaluate => integration_evaluate
     procedure :: make_iterations_list => integration_make_iterations_list
     procedure :: init_iteration_multipliers => &
          integration_init_iteration_multipliers
     procedure :: apply_call_multipliers => integration_apply_call_multipliers
     procedure :: init => integration_init
     procedure :: integrate => integration_integrate
     procedure :: integrate_dummy => integration_integrate_dummy
     procedure :: sampler_test => integration_sampler_test
     procedure :: get_process_ptr => integration_get_process_ptr
  end type integration_t


  interface
    module subroutine integration_init_process (intg, local)
      class(integration_t), intent(inout) :: intg
      type(rt_data_t), intent(inout), target :: local
    end subroutine integration_init_process
    module subroutine integration_setup_process &
         (intg, local, verbose, init_only)
      class(integration_t), intent(inout) :: intg
      type(rt_data_t), intent(inout), target :: local
      logical, intent(in), optional :: verbose
      logical, intent(in), optional :: init_only
    end subroutine integration_setup_process
    module subroutine integration_evaluate &
         (intg, process_instance, i_mci, pass, it_list, pacify)
      class(integration_t), intent(inout) :: intg
      type(process_instance_t), intent(inout), target :: process_instance
      integer, intent(in) :: i_mci
      integer, intent(in) :: pass
      type(iterations_list_t), intent(in) :: it_list
      logical, intent(in), optional :: pacify
    end subroutine integration_evaluate
    module subroutine integration_make_iterations_list (intg, it_list)
      class(integration_t), intent(in) :: intg
      type(iterations_list_t), intent(out) :: it_list
    end subroutine integration_make_iterations_list
    module subroutine integration_init_iteration_multipliers (intg, local)
      class(integration_t), intent(inout) :: intg
      type(rt_data_t), intent(in) :: local
    end subroutine integration_init_iteration_multipliers
    module subroutine integration_apply_call_multipliers &
         (intg, n_pass, i_component, it_list)
      class(integration_t), intent(in) :: intg
      integer, intent(in) :: n_pass, i_component
      type(iterations_list_t), intent(inout) :: it_list
    end subroutine integration_apply_call_multipliers
    module subroutine integration_init &
         (intg, process_id, local, global, local_stack, init_only)
      class(integration_t), intent(out) :: intg
      type(string_t), intent(in) :: process_id
      type(rt_data_t), intent(inout), target :: local
      type(rt_data_t), intent(inout), optional, target :: global
      logical, intent(in), optional :: init_only
      logical, intent(in), optional :: local_stack
    end subroutine integration_init
    module subroutine integration_integrate (intg, local, eff_reset)
      class(integration_t), intent(inout) :: intg
      type(rt_data_t), intent(in), target :: local
      logical, intent(in), optional :: eff_reset
    end subroutine integration_integrate
    module subroutine integration_integrate_dummy (intg)
      class(integration_t), intent(inout) :: intg
    end subroutine integration_integrate_dummy
    module subroutine integration_sampler_test (intg)
      class(integration_t), intent(inout) :: intg
    end subroutine integration_sampler_test
    module function integration_get_process_ptr (intg) result (ptr)
      class(integration_t), intent(in) :: intg
      type(process_t), pointer :: ptr
    end function integration_get_process_ptr
    module subroutine integrate_process &
         (process_id, local, global, local_stack, init_only, eff_reset)
      type(string_t), intent(in) :: process_id
      type(rt_data_t), intent(inout), target :: local
      type(rt_data_t), intent(inout), optional, target :: global
      logical, intent(in), optional :: local_stack, init_only, eff_reset
    end subroutine integrate_process
  end interface

contains

  subroutine integration_create_process (intg, process_id, global)
    class(integration_t), intent(out) :: intg
    type(rt_data_t), intent(inout), optional, target :: global
    type(string_t), intent(in) :: process_id
    type(process_entry_t), pointer :: process_entry
    if (debug_on) call msg_debug (D_CORE, "integration_create_process")
    intg%process_id = process_id
    if (present (global)) then
       allocate (process_entry)
       intg%process => process_entry%process_t
       call global%process_stack%push (process_entry)
    else
       allocate (process_t :: intg%process)
    end if
  end subroutine integration_create_process


end module integrations
