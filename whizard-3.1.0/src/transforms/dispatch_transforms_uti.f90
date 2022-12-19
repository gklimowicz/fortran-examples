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

module dispatch_transforms_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use format_utils, only: write_separator
  use variables
  use event_base, only: event_callback_t
  use models, only: model_t, model_list_t
  use models, only: syntax_model_file_init, syntax_model_file_final
  use resonances, only: resonance_history_set_t
  use beam_structures, only: beam_structure_t
  use eio_base, only: eio_t
  use os_interface, only: os_data_t
  use event_transforms, only: evt_t
  use dispatch_transforms

  implicit none
  private

  public :: dispatch_transforms_1
  public :: dispatch_transforms_2

contains

  subroutine dispatch_transforms_1 (u)
    integer, intent(in) :: u
    type(var_list_t) :: var_list
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(os_data_t) :: os_data
    class(event_callback_t), allocatable :: event_callback
    class(eio_t), allocatable :: eio

    write (u, "(A)")  "* Test output: dispatch_transforms_1"
    write (u, "(A)")  "*   Purpose: allocate an event I/O (eio) stream"
    write (u, "(A)")

    call var_list%init_defaults (0)
    call os_data%init ()
    call syntax_model_file_init ()
    call model_list%read_model (var_str ("SM_hadrons"), &
         var_str ("SM_hadrons.mdl"), os_data, model)

    write (u, "(A)")  "* Allocate as raw"
    write (u, "(A)")

    call dispatch_eio (eio, var_str ("raw"), var_list, &
         model, event_callback)

    call eio%write (u)

    call eio%final ()
    deallocate (eio)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate as checkpoints:"
    write (u, "(A)")

    call dispatch_eio (eio, var_str ("checkpoint"), var_list, &
         model, event_callback)

    call eio%write (u)

    call eio%final ()
    deallocate (eio)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate as LHEF:"
    write (u, "(A)")

    call var_list%set_string (var_str ("$lhef_extension"), &
         var_str ("lhe_custom"), is_known = .true.)
    call dispatch_eio (eio, var_str ("lhef"), var_list, &
         model, event_callback)

    call eio%write (u)

    call eio%final ()
    deallocate (eio)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate as HepMC:"
    write (u, "(A)")

    call dispatch_eio (eio, var_str ("hepmc"), var_list, &
         model, event_callback)

    call eio%write (u)

    call eio%final ()
    deallocate (eio)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate as weight_stream"
    write (u, "(A)")

    call dispatch_eio (eio, var_str ("weight_stream"), var_list, &
         model, event_callback)

    call eio%write (u)

    call eio%final ()
    deallocate (eio)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate as debug format"
    write (u, "(A)")

    call var_list%set_log (var_str ("?debug_verbose"), &
         .false., is_known = .true.)
    call dispatch_eio (eio, var_str ("debug"), var_list, &
         model, event_callback)

    call eio%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eio%final ()
    call var_list%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_transforms_1"

  end subroutine dispatch_transforms_1

  subroutine dispatch_transforms_2 (u)
    integer, intent(in) :: u
    type(var_list_t), target :: var_list
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(os_data_t) :: os_data
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(beam_structure_t) :: beam_structure
    class(evt_t), pointer :: evt

    write (u, "(A)")  "* Test output: dispatch_transforms_2"
    write (u, "(A)")  "*   Purpose: configure event transform"
    write (u, "(A)")

    call syntax_model_file_init ()
    call var_list%init_defaults (0)
    call os_data%init ()
    call model_list%read_model (var_str ("SM_hadrons"), &
         var_str ("SM_hadrons.mdl"), os_data, model)

    write (u, "(A)")  "* Resonance insertion"
    write (u, "(A)")

    call var_list%set_log (var_str ("?resonance_history"), .true., &
         is_known = .true.)
    call dispatch_evt_resonance (evt, var_list, &
         res_history_set, &
         var_str ("foo_R"))
    call evt%write (u, verbose = .true., more_verbose = .true.)

    call evt%final ()
    deallocate (evt)

    write (u, "(A)")
    write (u, "(A)")  "* ISR handler"
    write (u, "(A)")

    call var_list%set_log (var_str ("?isr_handler"), .true., &
         is_known = .true.)
    call var_list%set_log (var_str ("?epa_handler"), .false., &
         is_known = .true.)
    call var_list%set_string (var_str ("$isr_handler_mode"), &
         var_str ("recoil"), &
         is_known = .true.)
    call var_list%set_real (var_str ("sqrts"), 100._default, &
         is_known = .true.)
    call var_list%set_real (var_str ("isr_mass"), 511.e-6_default, &
         is_known = .true.)
    call dispatch_evt_isr_epa_handler (evt, var_list)
    call evt%write (u, verbose = .true., more_verbose = .true.)

    call evt%final ()
    deallocate (evt)

    write (u, "(A)")
    write (u, "(A)")  "* EPA handler"
    write (u, "(A)")

    call var_list%set_log (var_str ("?isr_handler"), .false., &
         is_known = .true.)
    call var_list%set_log (var_str ("?epa_handler"), .true., &
         is_known = .true.)
    call var_list%set_string (var_str ("$epa_handler_mode"), &
         var_str ("recoil"), &
         is_known = .true.)
    call var_list%set_real (var_str ("sqrts"), 100._default, &
         is_known = .true.)
    call var_list%set_real (var_str ("epa_mass"), 511.e-6_default, &
         is_known = .true.)
    call dispatch_evt_isr_epa_handler (evt, var_list)
    call evt%write (u, verbose = .true., more_verbose = .true.)

    call evt%final ()
    deallocate (evt)

    write (u, "(A)")
    write (u, "(A)")  "* Partonic decays"
    write (u, "(A)")

    call dispatch_evt_decay (evt, var_list)
    call evt%write (u, verbose = .true., more_verbose = .true.)

    call evt%final ()
    deallocate (evt)

    write (u, "(A)")
    write (u, "(A)")  "* Shower"
    write (u, "(A)")

    call var_list%set_log (var_str ("?allow_shower"), .true., &
         is_known = .true.)
    call var_list%set_string (var_str ("$shower_method"), &
         var_str ("WHIZARD"), is_known = .true.)
    call dispatch_evt_shower (evt, var_list, model, &
         model, os_data, beam_structure)
    call evt%write (u)
    call write_separator (u, 2)

    call evt%final ()
    deallocate (evt)

    call var_list%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_transforms_2"

  end subroutine dispatch_transforms_2


end module dispatch_transforms_uti
