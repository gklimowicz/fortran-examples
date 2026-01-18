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

module shower

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use pdf

  use shower_base
  use matching_base

  use sm_qcd
  use model_data

  use event_transforms
  use models
  use process, only: process_t
  use instances, only: process_instance_t
  use process_stacks

  implicit none
  private

  public :: evt_shower_t

  logical, parameter :: POWHEG_TESTING = .true.


  type, extends (evt_t) :: evt_shower_t
     class(shower_base_t), allocatable :: shower
     class(matching_t), allocatable :: matching
     type(model_t), pointer :: model_hadrons => null ()
     type(qcd_t) :: qcd
     type(pdf_data_t) :: pdf_data
     type(os_data_t) :: os_data
     logical :: is_first_event
   contains
     procedure :: write_name => evt_shower_write_name
     procedure :: write => evt_shower_write
     procedure :: connect => evt_shower_connect
     procedure :: init => evt_shower_init
     procedure :: make_rng => evt_shower_make_rng
     procedure :: prepare_new_event => evt_shower_prepare_new_event
     procedure :: first_event => evt_shower_first_event
     procedure :: generate_weighted => evt_shower_generate_weighted
     procedure :: make_particle_set => evt_shower_make_particle_set
     procedure :: contains_powheg_matching => evt_shower_contains_powheg_matching
     procedure :: disable_powheg_matching => evt_shower_disable_powheg_matching
     procedure :: enable_powheg_matching => evt_shower_enable_powheg_matching
     procedure :: final => evt_shower_final
  end type evt_shower_t


  interface
    module subroutine evt_shower_write_name (evt, unit)
      class(evt_shower_t), intent(in) :: evt
      integer, intent(in), optional :: unit
    end subroutine evt_shower_write_name
    module subroutine evt_shower_write &
         (evt, unit, verbose, more_verbose, testflag)
      class(evt_shower_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, more_verbose, testflag
    end subroutine evt_shower_write
    module subroutine evt_shower_connect &
         (evt, process_instance, model, process_stack)
      class(evt_shower_t), intent(inout), target :: evt
      type(process_instance_t), intent(in), target :: process_instance
      class(model_data_t), intent(in), target :: model
      type(process_stack_t), intent(in), optional :: process_stack
    end subroutine evt_shower_connect
    module subroutine evt_shower_init (evt, model_hadrons, os_data)
      class(evt_shower_t), intent(out) :: evt
      type(model_t), intent(in), target :: model_hadrons
      type(os_data_t), intent(in) :: os_data
    end subroutine evt_shower_init
    module subroutine evt_shower_make_rng (evt, process)
      class(evt_shower_t), intent(inout) :: evt
      type(process_t), intent(inout) :: process
    end subroutine evt_shower_make_rng
    module subroutine evt_shower_prepare_new_event (evt, i_mci, i_term)
      class(evt_shower_t), intent(inout) :: evt
      integer, intent(in) :: i_mci, i_term
    end subroutine evt_shower_prepare_new_event
    module subroutine evt_shower_first_event (evt)
      class(evt_shower_t), intent(inout) :: evt
    end subroutine evt_shower_first_event
    module subroutine evt_shower_generate_weighted (evt, probability)
      class(evt_shower_t), intent(inout) :: evt
      real(default), intent(inout) :: probability
    end subroutine evt_shower_generate_weighted
    module subroutine evt_shower_make_particle_set &
         (evt, factorization_mode, keep_correlations, r)
      class(evt_shower_t), intent(inout) :: evt
      integer, intent(in) :: factorization_mode
      logical, intent(in) :: keep_correlations
      real(default), dimension(:), intent(in), optional :: r
    end subroutine evt_shower_make_particle_set
    module function evt_shower_contains_powheg_matching (evt) result (val)
       logical :: val
       class(evt_shower_t), intent(in) :: evt
    end function evt_shower_contains_powheg_matching
    module subroutine evt_shower_disable_powheg_matching (evt)
       class(evt_shower_t), intent(inout) :: evt
    end subroutine evt_shower_disable_powheg_matching
    module subroutine evt_shower_enable_powheg_matching (evt)
       class(evt_shower_t), intent(inout) :: evt
    end subroutine evt_shower_enable_powheg_matching
    module subroutine evt_shower_final (evt)
      class(evt_shower_t), intent(inout) :: evt
    end subroutine evt_shower_final
  end interface

end module shower
