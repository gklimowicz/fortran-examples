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

module event_transforms_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use format_utils, only: write_separator
  use os_interface
  use sm_qcd
  use models
  use state_matrices, only: FM_IGNORE_HELICITY
  use interactions, only: reset_interaction_counter
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

  use event_transforms

  use rng_base_ut, only: rng_test_factory_t

  implicit none
  private

  public :: event_transforms_1

contains

  subroutine event_transforms_1 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(model_t), target :: model
    type(process_library_t), target :: lib
    type(string_t) :: libname, procname1
    class(phs_config_t), allocatable :: phs_config_template
    real(default) :: sqrts
    type(process_t), allocatable, target :: process
    type(process_instance_t), allocatable, target :: process_instance
    class(evt_t), allocatable :: evt
    integer :: factorization_mode
    logical :: keep_correlations

    write (u, "(A)")  "* Test output: event_transforms_1"
    write (u, "(A)")  "*   Purpose: handle trivial transform"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment and parent process"
    write (u, "(A)")

    call os_data%init ()

    libname = "event_transforms_1_lib"
    procname1 = "event_transforms_1_p"

    call prc_test_create_library (libname, lib, &
         scattering = .true., procname1 = procname1)
    call reset_interaction_counter ()

    call model%init_test ()

    allocate (process)
    call process%init (procname1, lib, os_data, model)
    call process%setup_test_cores ()

    allocate (phs_single_config_t :: phs_config_template)
    call process%init_components (phs_config_template)

    sqrts = 1000
    call process%setup_beams_sqrts (sqrts, i_core = 1)
    call process%configure_phs ()
    call process%setup_mci (dispatch_mci_test_midpoint)
    call process%setup_terms ()

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%integrate (1, n_it=1, n_calls=100)
    call process%final_integration (1)
    call process_instance%final ()
    deallocate (process_instance)

    allocate (process_instance)
    call process_instance%init (process)
    call process_instance%setup_event_data ()
    call process_instance%init_simulation (1)

    write (u, "(A)")  "* Initialize trivial event transform"
    write (u, "(A)")

    allocate (evt_trivial_t :: evt)
    call evt%connect (process_instance, process%get_model_ptr ())

    write (u, "(A)")  "* Generate event and subsequent transform"
    write (u, "(A)")

    call process_instance%generate_unweighted_event (1)
    call process_instance%evaluate_event_data ()

    call evt%prepare_new_event (1, 1)
    call evt%generate_unweighted ()

    call write_separator (u, 2)
    call evt%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Obtain particle set"
    write (u, "(A)")

    factorization_mode = FM_IGNORE_HELICITY
    keep_correlations = .false.

    call evt%make_particle_set (factorization_mode, keep_correlations)

    call write_separator (u, 2)
    call evt%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt%final ()
    call process_instance%final ()
    deallocate (process_instance)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: event_transforms_1"

  end subroutine event_transforms_1


  subroutine dispatch_mci_test_midpoint (mci, var_list, process_id, is_nlo)
    use variables, only: var_list_t
    class(mci_t), allocatable, intent(out) :: mci
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_id
    logical, intent(in), optional :: is_nlo
    allocate (mci_midpoint_t :: mci)
  end subroutine dispatch_mci_test_midpoint


end module event_transforms_uti

