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

module hep_events

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use polarizations
  use model_data
  use particles
  use hepmc_interface
  use lcio_interface
  use event_base

  implicit none
  private

  public :: hepeup_from_event
  public :: hepeup_to_event
  public :: hepevt_from_event
  public :: hepmc_event_from_particle_set
  public :: hepmc_event_to_particle_set
  public :: hepmc_to_event
  public :: particle_to_lcio
  public :: particle_from_lcio_particle
  public :: lcio_event_from_particle_set
  public :: lcio_event_to_particle_set
  public :: lcio_to_event

  interface
    module subroutine hepeup_from_event &
         (event, keep_beams, keep_remnants, process_index)
      class(generic_event_t), intent(in), target :: event
      logical, intent(in), optional :: keep_beams
      logical, intent(in), optional :: keep_remnants
      integer, intent(in), optional :: process_index
    end subroutine hepeup_from_event
    module subroutine hepeup_to_event &
         (event, fallback_model, process_index, recover_beams, &
         use_alpha_s, use_scale)
      class(generic_event_t), intent(inout), target :: event
      class(model_data_t), intent(in), target :: fallback_model
      integer, intent(out), optional :: process_index
      logical, intent(in), optional :: recover_beams
      logical, intent(in), optional :: use_alpha_s
      logical, intent(in), optional :: use_scale
    end subroutine hepeup_to_event
    module subroutine hepevt_from_event  &
           (event, process_index, i_evt, keep_beams, keep_remnants, &
            ensure_order, fill_hepev4)
      class(generic_event_t), intent(in), target :: event
      integer, intent(in), optional :: i_evt, process_index
      logical, intent(in), optional :: keep_beams
      logical, intent(in), optional :: keep_remnants
      logical, intent(in), optional :: ensure_order
      logical, intent(in), optional :: fill_hepev4
    end subroutine hepevt_from_event
    module subroutine hepmc_event_from_particle_set &
           (evt, particle_set, cross_section, error, color)
      type(hepmc_event_t), intent(inout) :: evt
      type(particle_set_t), intent(in) :: particle_set
      real(default), intent(in), optional :: cross_section, error
      logical, intent(in), optional :: color
    end subroutine hepmc_event_from_particle_set
    module subroutine hepmc_event_to_particle_set &
         (particle_set, evt, model, fallback_model, polarization)
      type(particle_set_t), intent(inout), target :: particle_set
      type(hepmc_event_t), intent(in) :: evt
      class(model_data_t), intent(in), target :: model, fallback_model
      integer, intent(in) :: polarization
    end subroutine hepmc_event_to_particle_set
    module subroutine hepmc_to_event &
         (event, hepmc_event, fallback_model, process_index, &
         recover_beams, use_alpha_s, use_scale)
      class(generic_event_t), intent(inout), target :: event
      type(hepmc_event_t), intent(inout) :: hepmc_event
      class(model_data_t), intent(in), target :: fallback_model
      integer, intent(out), optional :: process_index
      logical, intent(in), optional :: recover_beams
      logical, intent(in), optional :: use_alpha_s
      logical, intent(in), optional :: use_scale
    end subroutine hepmc_to_event
    module subroutine particle_to_lcio (prt, lprt)
      type(particle_t), intent(in) :: prt
      type(lcio_particle_t), intent(out) :: lprt
    end subroutine particle_to_lcio
    module subroutine particle_from_lcio_particle &
       (prt, lprt, model, fallback_model, daughters, parents, polarization)
      type(particle_t), intent(out) :: prt
      type(lcio_particle_t), intent(in) :: lprt
      type(model_data_t), intent(in), target :: model
      type(model_data_t), intent(in), target :: fallback_model
      integer, dimension(:), intent(in) :: daughters, parents
      integer, intent(in) :: polarization
    end subroutine particle_from_lcio_particle
    module subroutine lcio_event_from_particle_set (evt, particle_set)
      type(lcio_event_t), intent(inout) :: evt
      type(particle_set_t), intent(in) :: particle_set
    end subroutine lcio_event_from_particle_set
    module subroutine lcio_event_to_particle_set &
         (particle_set, evt, model, fallback_model, polarization)
      type(particle_set_t), intent(inout), target :: particle_set
      type(lcio_event_t), intent(in) :: evt
      class(model_data_t), intent(in), target :: model, fallback_model
      integer, intent(in) :: polarization
    end subroutine lcio_event_to_particle_set
    module subroutine lcio_to_event &
         (event, lcio_event, fallback_model, process_index, recover_beams, &
         use_alpha_s, use_scale)
      class(generic_event_t), intent(inout), target :: event
      type(lcio_event_t), intent(inout) :: lcio_event
      class(model_data_t), intent(in), target :: fallback_model
      integer, intent(out), optional :: process_index
      logical, intent(in), optional :: recover_beams
      logical, intent(in), optional :: use_alpha_s
      logical, intent(in), optional :: use_scale
    end subroutine lcio_to_event
  end interface

end module hep_events
