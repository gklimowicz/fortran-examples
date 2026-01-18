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

module event_transforms

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use quantum_numbers, only: quantum_numbers_t
  use interactions
  use particles
  use model_data
  use rng_base
  use process, only: process_t
  use instances, only: process_instance_t
  use process_stacks

  implicit none
  private

  public :: evt_t
  public :: make_factorized_particle_set
  public :: evt_trivial_t

  type, abstract :: evt_t
     type(process_t), pointer :: process => null ()
     type(process_instance_t), pointer :: process_instance => null ()
     class(model_data_t), pointer :: model => null ()
     class(rng_t), allocatable :: rng
     integer :: rejection_count = 0
     logical :: particle_set_exists = .false.
     type(particle_set_t) :: particle_set
     class(evt_t), pointer :: previous => null ()
     class(evt_t), pointer :: next => null ()
     logical :: only_weighted_events = .false.
   contains
     procedure :: final => evt_final
     procedure :: base_final => evt_final
     procedure (evt_write_name), deferred :: write_name
     procedure (evt_write), deferred :: write
     procedure :: base_write => evt_base_write
     procedure :: connect => evt_connect
     procedure :: base_connect => evt_connect
     procedure :: reset => evt_reset
     procedure :: base_reset => evt_reset
     procedure (evt_prepare_new_event), deferred :: prepare_new_event
     procedure (evt_generate_weighted), deferred :: generate_weighted
     procedure :: generate_unweighted => evt_generate_unweighted
     procedure :: base_generate_unweighted => evt_generate_unweighted
     procedure (evt_make_particle_set), deferred :: make_particle_set
     procedure :: set_particle_set => evt_set_particle_set
     procedure :: factorize_interactions => evt_factorize_interactions
     procedure :: tag_incoming => evt_tag_incoming
  end type evt_t

  type, extends (evt_t) :: evt_trivial_t
   contains
     procedure :: write_name => evt_trivial_write_name
     procedure :: write => evt_trivial_write
     procedure :: prepare_new_event => evt_trivial_prepare_new_event
     procedure :: generate_weighted => evt_trivial_generate_weighted
     procedure :: make_particle_set => evt_trivial_make_particle_set
  end type evt_trivial_t


  abstract interface
     subroutine evt_write_name (evt, unit)
       import
       class(evt_t), intent(in) :: evt
       integer, intent(in), optional :: unit
     end subroutine evt_write_name
  end interface

  abstract interface
     subroutine evt_write (evt, unit, verbose, more_verbose, testflag)
       import
       class(evt_t), intent(in) :: evt
       integer, intent(in), optional :: unit
       logical, intent(in), optional :: verbose, more_verbose, testflag
     end subroutine evt_write
  end interface

  interface
     subroutine evt_prepare_new_event (evt, i_mci, i_term)
       import
       class(evt_t), intent(inout) :: evt
       integer, intent(in) :: i_mci, i_term
     end subroutine evt_prepare_new_event
  end interface

  abstract interface
     subroutine evt_generate_weighted (evt, probability)
       import
       class(evt_t), intent(inout) :: evt
       real(default), intent(inout) :: probability
     end subroutine evt_generate_weighted
  end interface

  interface
     subroutine evt_make_particle_set &
          (evt, factorization_mode, keep_correlations, r)
       import
       class(evt_t), intent(inout) :: evt
       integer, intent(in) :: factorization_mode
       logical, intent(in) :: keep_correlations
       real(default), dimension(:), intent(in), optional :: r
     end subroutine evt_make_particle_set
  end interface


  interface
    module subroutine evt_final (evt)
      class(evt_t), intent(inout) :: evt
    end subroutine evt_final
    module subroutine evt_base_write (evt, unit, testflag, show_set)
      class(evt_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag, show_set
    end subroutine evt_base_write
    module subroutine evt_connect (evt, process_instance, model, process_stack)
      class(evt_t), intent(inout), target :: evt
      type(process_instance_t), intent(in), target :: process_instance
      class(model_data_t), intent(in), target :: model
      type(process_stack_t), intent(in), optional :: process_stack
    end subroutine evt_connect
    module subroutine evt_reset (evt)
      class(evt_t), intent(inout) :: evt
    end subroutine evt_reset
    module subroutine evt_generate_unweighted (evt)
      class(evt_t), intent(inout) :: evt
    end subroutine evt_generate_unweighted
    module subroutine evt_set_particle_set (evt, particle_set, i_mci, i_term)
      class(evt_t), intent(inout) :: evt
      type(particle_set_t), intent(in) :: particle_set
      integer, intent(in) :: i_term, i_mci
    end subroutine evt_set_particle_set
    module subroutine evt_factorize_interactions &
         (evt, int_matrix, int_flows, factorization_mode, &
         keep_correlations, r, qn_select)
      class(evt_t), intent(inout) :: evt
      type(interaction_t), intent(in), target :: int_matrix, int_flows
      integer, intent(in) :: factorization_mode
      logical, intent(in) :: keep_correlations
      real(default), dimension(:), intent(in), optional :: r
      type(quantum_numbers_t), dimension(:), intent(in), optional :: qn_select
    end subroutine evt_factorize_interactions
    module subroutine make_factorized_particle_set (evt, factorization_mode, &
           keep_correlations, r, ii_term, qn_select)
      class(evt_t), intent(inout) :: evt
      integer, intent(in) :: factorization_mode
      logical, intent(in) :: keep_correlations
      real(default), dimension(:), intent(in), optional :: r
      integer, intent(in), optional :: ii_term
      type(quantum_numbers_t), dimension(:), intent(in), optional :: qn_select
    end subroutine make_factorized_particle_set
    module subroutine evt_tag_incoming (evt)
      class(evt_t), intent(inout) :: evt
    end subroutine evt_tag_incoming
    module subroutine evt_trivial_write_name (evt, unit)
      class(evt_trivial_t), intent(in) :: evt
      integer, intent(in), optional :: unit
    end subroutine evt_trivial_write_name
    module subroutine evt_trivial_write &
         (evt, unit, verbose, more_verbose, testflag)
      class(evt_trivial_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, more_verbose, testflag
    end subroutine evt_trivial_write
    module subroutine evt_trivial_prepare_new_event (evt, i_mci, i_term)
      class(evt_trivial_t), intent(inout) :: evt
      integer, intent(in) :: i_mci, i_term
    end subroutine evt_trivial_prepare_new_event
    module subroutine evt_trivial_generate_weighted (evt, probability)
      class(evt_trivial_t), intent(inout) :: evt
      real(default), intent(inout) :: probability
    end subroutine evt_trivial_generate_weighted
    module subroutine evt_trivial_make_particle_set &
         (evt, factorization_mode, keep_correlations, r)
      class(evt_trivial_t), intent(inout) :: evt
      integer, intent(in) :: factorization_mode
      logical, intent(in) :: keep_correlations
      real(default), dimension(:), intent(in), optional :: r
    end subroutine evt_trivial_make_particle_set
  end interface

end module event_transforms

