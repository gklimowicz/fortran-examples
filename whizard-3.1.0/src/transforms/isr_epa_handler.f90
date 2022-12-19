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

module isr_epa_handler

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz, only: vector4_t
  use lorentz, only: energy
  use lorentz, only: lorentz_transformation_t
  use lorentz, only: identity
  use lorentz, only: inverse
  use lorentz, only: operator(*)
  use flavors, only: flavor_t
  use particles, only: particle_t
  use rng_base, only: rng_t
  use event_transforms

  implicit none
  private

  public :: evt_isr_epa_t

  integer, parameter, public :: BEAM_RAD_NONE = 0
  integer, parameter, public :: BEAM_RAD_ISR = 1
  integer, parameter, public :: BEAM_RAD_EPA = 2

  integer, parameter, public :: ISR_TRIVIAL_COLLINEAR = 0
  integer, parameter, public :: ISR_PAIR_RECOIL = 1


  type, extends (evt_t) :: evt_isr_epa_t
     private
     integer :: mode = ISR_TRIVIAL_COLLINEAR
     logical :: isr_active = .false.
     logical :: epa_active = .false.
     real(default) :: isr_q_max = 0
     real(default) :: epa_q_max = 0
     real(default) :: isr_mass = 0
     real(default) :: epa_mass = 0
     logical :: isr_keep_mass = .true.
     real(default) :: sqrts = 0
     integer, dimension(2) :: rad_mode = BEAM_RAD_NONE
     real(default), dimension(2) :: q_max = 0
     real(default), dimension(2) :: m = 0
     real(default), dimension(2) :: xc = 0
     real(default), dimension(2) :: xcb = 0
     type(lorentz_transformation_t) :: lti = identity
     type(lorentz_transformation_t) :: lto = identity
     type(lorentz_transformation_t) :: lt  = identity
     integer, dimension(2) :: i_beam = 0
     type(particle_t), dimension(2) :: beam
     type(vector4_t), dimension(2) :: pi
     integer, dimension(2) :: i_radiated = 0
     type(particle_t), dimension(2) :: radiated
     type(vector4_t), dimension(2) :: ki
     type(vector4_t), dimension(2) :: km
     integer, dimension(2) :: i_parton = 0
     type(particle_t), dimension(2) :: parton
     type(vector4_t), dimension(2) :: qi
     type(vector4_t), dimension(2) :: qm
     type(vector4_t), dimension(2) :: qo
   contains
     procedure :: get_mode_string => evt_isr_epa_get_mode_string
     procedure :: set_mode_string => evt_isr_epa_set_mode_string
     procedure :: write_name => evt_isr_epa_write_name
     procedure :: write_mode => evt_isr_epa_write_mode
     procedure :: write_input => evt_isr_epa_write_input
     procedure :: write_data => evt_isr_epa_write_data
     procedure :: write => evt_isr_epa_write
     procedure :: import_rng => evt_isr_epa_import_rng
     procedure :: set_data_isr => evt_isr_epa_set_data_isr
     procedure :: set_data_epa => evt_isr_epa_set_data_epa
     procedure, private :: identify_radiated
     procedure, private :: identify_partons
     procedure :: check_radiation => evt_isr_epa_check_radiation
     procedure :: set_recoil_parameters => evt_isr_epa_set_recoil_parameters
     procedure, private :: boost_to_cm
     procedure, private :: infer_x
     procedure, private :: generate_recoil => evt_generate_recoil
     procedure, private :: replace_radiated
     procedure, private :: replace_partons
     procedure :: transform_outgoing => evt_transform_outgoing
     procedure :: generate_weighted => &
          evt_isr_epa_generate_weighted
     procedure :: make_particle_set => &
          evt_isr_epa_make_particle_set
     procedure :: prepare_new_event => &
          evt_isr_epa_prepare_new_event
  end type evt_isr_epa_t


  interface
    module function evt_isr_epa_get_mode_string (evt) result (string)
      type(string_t) :: string
      class(evt_isr_epa_t), intent(in) :: evt
    end function evt_isr_epa_get_mode_string
    module subroutine evt_isr_epa_set_mode_string (evt, string)
      class(evt_isr_epa_t), intent(inout) :: evt
      type(string_t), intent(in) :: string
    end subroutine evt_isr_epa_set_mode_string
    module subroutine evt_isr_epa_write_name (evt, unit)
      class(evt_isr_epa_t), intent(in) :: evt
      integer, intent(in), optional :: unit
    end subroutine evt_isr_epa_write_name
    module subroutine evt_isr_epa_write_mode (evt, unit)
      class(evt_isr_epa_t), intent(in) :: evt
      integer, intent(in), optional :: unit
    end subroutine evt_isr_epa_write_mode
    module subroutine evt_isr_epa_write_input (evt, unit, testflag)
      class(evt_isr_epa_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine evt_isr_epa_write_input
    module subroutine evt_isr_epa_write_data (evt, unit, testflag)
      class(evt_isr_epa_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine evt_isr_epa_write_data
    module subroutine evt_isr_epa_write &
         (evt, unit, verbose, more_verbose, testflag)
      class(evt_isr_epa_t), intent(in) :: evt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose, more_verbose, testflag
    end subroutine evt_isr_epa_write
    module subroutine evt_isr_epa_import_rng (evt, rng)
      class(evt_isr_epa_t), intent(inout) :: evt
      class(rng_t), allocatable, intent(inout) :: rng
    end subroutine evt_isr_epa_import_rng
    module subroutine evt_isr_epa_set_data_isr (evt, sqrts, q_max, m, keep_mass)
      class(evt_isr_epa_t), intent(inout) :: evt
      real(default), intent(in) :: sqrts
      real(default), intent(in) :: q_max
      real(default), intent(in) :: m
      logical, intent(in) :: keep_mass
    end subroutine evt_isr_epa_set_data_isr
    module subroutine evt_isr_epa_set_data_epa (evt, sqrts, q_max, m)
      class(evt_isr_epa_t), intent(inout) :: evt
      real(default), intent(in) :: sqrts
      real(default), intent(in) :: q_max
      real(default), intent(in) :: m
    end subroutine evt_isr_epa_set_data_epa
    module subroutine identify_radiated (evt)
      class(evt_isr_epa_t), intent(inout) :: evt
    end subroutine identify_radiated
    module subroutine identify_partons (evt)
      class(evt_isr_epa_t), intent(inout) :: evt
    end subroutine identify_partons
    module subroutine evt_isr_epa_check_radiation (evt)
      class(evt_isr_epa_t), intent(inout) :: evt
    end subroutine evt_isr_epa_check_radiation
    module subroutine evt_isr_epa_set_recoil_parameters (evt)
      class(evt_isr_epa_t), intent(inout) :: evt
    end subroutine evt_isr_epa_set_recoil_parameters
    module subroutine boost_to_cm (evt)
      class(evt_isr_epa_t), intent(inout) :: evt
    end subroutine boost_to_cm
    module subroutine infer_x (evt)
      class(evt_isr_epa_t), intent(inout) :: evt
    end subroutine infer_x
    module subroutine evt_generate_recoil (evt, ok)
      class(evt_isr_epa_t), intent(inout) :: evt
      logical, intent(out) :: ok
    end subroutine evt_generate_recoil
    module subroutine replace_radiated (evt)
      class(evt_isr_epa_t), intent(inout) :: evt
    end subroutine replace_radiated
    module subroutine replace_partons (evt)
      class(evt_isr_epa_t), intent(inout) :: evt
    end subroutine replace_partons
    module subroutine evt_transform_outgoing (evt)
      class(evt_isr_epa_t), intent(inout) :: evt
    end subroutine evt_transform_outgoing
    module subroutine evt_isr_epa_generate_weighted (evt, probability)
      class(evt_isr_epa_t), intent(inout) :: evt
      real(default), intent(inout) :: probability
    end subroutine evt_isr_epa_generate_weighted
    module subroutine evt_isr_epa_make_particle_set &
         (evt, factorization_mode, keep_correlations, r)
      class(evt_isr_epa_t), intent(inout) :: evt
      integer, intent(in) :: factorization_mode
      logical, intent(in) :: keep_correlations
      real(default), dimension(:), intent(in), optional :: r
    end subroutine evt_isr_epa_make_particle_set
    module subroutine evt_isr_epa_prepare_new_event (evt, i_mci, i_term)
      class(evt_isr_epa_t), intent(inout) :: evt
      integer, intent(in) :: i_mci, i_term
    end subroutine evt_isr_epa_prepare_new_event
  end interface

end module isr_epa_handler
