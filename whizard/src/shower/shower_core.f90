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

module shower_core

  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use os_interface
  use lorentz
  use particles
  use model_data
  use pdf
  use rng_base
  use shower_base
  use shower_partons
  use muli, only: muli_t
  use tauola_interface

  implicit none
  private

  public :: shower_interaction_t
  public :: shower_t
  public :: shower_interaction_get_s

  type :: shower_interaction_t
     type(parton_pointer_t), dimension(:), allocatable :: partons
  end type shower_interaction_t

  type :: shower_interaction_pointer_t
     type(shower_interaction_t), pointer :: i => null ()
  end type shower_interaction_pointer_t

  type, extends (shower_base_t) :: shower_t
     type(shower_interaction_pointer_t), dimension(:), allocatable :: &
          interactions
     type(parton_pointer_t), dimension(:), allocatable :: partons
     type(muli_t) :: mi
     integer :: next_free_nr
     integer :: next_color_nr
     logical :: valid
   contains
     procedure :: init => shower_init
     procedure :: prepare_new_event => shower_prepare_new_event
     procedure :: activate_multiple_interactions => &
          shower_activate_multiple_interactions
     procedure :: import_particle_set => shower_import_particle_set
     procedure :: generate_emissions => shower_generate_emissions
     procedure :: make_particle_set => shower_make_particle_set
     procedure :: add_interaction_2ton => shower_add_interaction_2ton
     procedure :: simulate_no_isr_shower => shower_simulate_no_isr_shower
     procedure :: simulate_no_fsr_shower => shower_simulate_no_fsr_shower
     procedure :: sort_partons => shower_sort_partons
     procedure :: cleanup => shower_cleanup
     procedure :: get_next_free_nr => shower_get_next_free_nr
     procedure :: update_max_color_nr => shower_update_max_color_nr
     procedure :: get_next_color_nr => shower_get_next_color_nr
     procedure :: add_child => shower_add_child
     procedure :: add_parent => shower_add_parent
     procedure :: get_nr_of_partons => shower_get_nr_of_partons
     procedure :: get_final_colored_ME_momenta => &
          shower_get_final_colored_ME_momenta
     procedure :: update_beamremnants => shower_update_beamremnants
     procedure :: boost_to_CMframe => shower_boost_to_CMframe
     procedure :: boost_to_labframe => shower_boost_to_labframe
     procedure :: rotate_to_z => shower_rotate_to_z
     procedure :: generate_primordial_kt => shower_generate_primordial_kt
     procedure :: write => shower_write
     procedure :: combine_with_particle_set => shower_combine_with_particle_set
     procedure :: write_lhef => shower_write_lhef
     procedure :: generate_next_isr_branching_veto => &
          shower_generate_next_isr_branching_veto
     procedure :: find_recoiler => shower_find_recoiler
     procedure :: generate_next_isr_branching => &
          shower_generate_next_isr_branching
     procedure :: generate_fsr_for_isr_partons => &
          shower_generate_fsr_for_partons_emitted_in_ISR
     procedure :: execute_next_isr_branching => shower_execute_next_isr_branching
     procedure :: get_ISR_scale => shower_get_ISR_scale
     procedure :: set_max_isr_scale => shower_set_max_isr_scale
     procedure :: interaction_generate_fsr_2ton => &
          shower_interaction_generate_fsr_2ton
     procedure :: parton_pointer_array_generate_fsr_recursive => &
          shower_parton_pointer_array_generate_fsr_recursive
     procedure :: parton_pointer_array_generate_fsr => &
          shower_parton_pointer_array_generate_fsr
     procedure :: get_pdf => shower_get_pdf
     procedure :: get_xpdf => shower_get_xpdf
     procedure :: converttopythia => shower_converttopythia
  end type shower_t




  interface
    module subroutine shower_init &
         (shower, settings, taudec_settings, pdf_data, os_data)
      class(shower_t), intent(out) :: shower
      type(shower_settings_t), intent(in) :: settings
      type(taudec_settings_t), intent(in) :: taudec_settings
      type(pdf_data_t), intent(in) :: pdf_data
      type(os_data_t), intent(in) :: os_data
    end subroutine shower_init
    module subroutine shower_prepare_new_event (shower, fac_scale, alpha_s)
      class(shower_t), intent(inout) :: shower
      real(default), intent(in) :: fac_scale, alpha_s
    end subroutine shower_prepare_new_event
    module subroutine shower_activate_multiple_interactions (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_activate_multiple_interactions
    module subroutine shower_import_particle_set (shower, particle_set)
      class(shower_t), target, intent(inout) :: shower
      type(particle_set_t), intent(in) :: particle_set
    end subroutine shower_import_particle_set
    module subroutine shower_generate_emissions &
         (shower, valid, number_of_emissions)
      class(shower_t), intent(inout), target :: shower
      logical, intent(out) :: valid
      integer, optional, intent(in) :: number_of_emissions
    end subroutine shower_generate_emissions
    module subroutine shower_make_particle_set &
         (shower, particle_set, model, model_hadrons)
      class(shower_t), intent(in) :: shower
      type(particle_set_t), intent(inout) :: particle_set
      class(model_data_t), intent(in), target :: model
      class(model_data_t), intent(in), target :: model_hadrons
    end subroutine shower_make_particle_set
    module subroutine shower_add_interaction_2ton (shower, partons)
      class(shower_t), intent(inout) :: shower
      type(parton_pointer_t), intent(in), dimension(:), allocatable :: partons
    end subroutine shower_add_interaction_2ton
    module subroutine shower_simulate_no_isr_shower (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_simulate_no_isr_shower
    module subroutine shower_simulate_no_fsr_shower (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_simulate_no_fsr_shower
    module subroutine shower_sort_partons (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_sort_partons
    module subroutine shower_cleanup (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_cleanup
    module function shower_get_next_free_nr (shower) result (next_number)
      class(shower_t), intent(inout) :: shower
      integer :: next_number
    end function shower_get_next_free_nr
    pure module subroutine shower_update_max_color_nr (shower, index)
      class(shower_t), intent(inout) :: shower
      integer, intent(in) :: index
    end subroutine shower_update_max_color_nr
    module function shower_get_next_color_nr (shower) result (next_color)
      class(shower_t), intent(inout) :: shower
      integer :: next_color
    end function shower_get_next_color_nr
    module subroutine shower_add_child (shower, prt, child)
      class(shower_t), intent(inout) :: shower
      type(parton_t), pointer :: prt
      integer, intent(in) :: child
    end subroutine shower_add_child
    module subroutine shower_add_parent (shower, prt)
      class(shower_t), intent(inout) :: shower
      type(parton_t), intent(inout), target :: prt
    end subroutine shower_add_parent
    module function shower_get_nr_of_partons (shower, mine, &
         include_remnants, no_hard_prts, only_colored) result (nr)
      class(shower_t), intent(in) :: shower
      real(default), intent(in), optional :: mine
      logical, intent(in), optional :: include_remnants, no_hard_prts, &
           only_colored
      integer :: nr
    end function shower_get_nr_of_partons
    module subroutine shower_get_final_colored_ME_momenta (shower, momenta)
      class(shower_t), intent(in) :: shower
      type(vector4_t), dimension(:), allocatable, intent(out) :: momenta
    end subroutine shower_get_final_colored_ME_momenta
    module function shower_interaction_get_s (interaction) result (s)
      type(shower_interaction_t), intent(in) :: interaction
      real(default) :: s
    end function shower_interaction_get_s
    module subroutine shower_update_beamremnants (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_update_beamremnants
    module subroutine shower_boost_to_CMframe (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_boost_to_CMframe
    module subroutine shower_boost_to_labframe (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_boost_to_labframe
    module subroutine shower_rotate_to_z (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_rotate_to_z
    module subroutine shower_generate_primordial_kt (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_generate_primordial_kt
    module subroutine shower_write (shower, unit)
      class(shower_t), intent(in) :: shower
      integer, intent(in), optional :: unit
    end subroutine shower_write
    module subroutine shower_combine_with_particle_set (shower, particle_set, &
         model_in, model_hadrons)
      class(shower_t), intent(in) :: shower
      type(particle_set_t), intent(inout) :: particle_set
      class(model_data_t), intent(in), target :: model_in
      class(model_data_t), intent(in), target :: model_hadrons
    end subroutine shower_combine_with_particle_set
    module subroutine shower_write_lhef (shower, unit)
      class(shower_t), intent(in) :: shower
      integer, intent(in), optional :: unit
    end subroutine shower_write_lhef
    module function shower_generate_next_isr_branching_veto &
         (shower) result (next_brancher)
      class(shower_t), intent(inout) :: shower
      type(parton_pointer_t) :: next_brancher
    end function shower_generate_next_isr_branching_veto
    module function shower_find_recoiler (shower, prt) result(recoiler)
      class(shower_t), intent(inout) :: shower
      type(parton_t), intent(inout), target :: prt
      type(parton_t), pointer :: recoiler
    end function shower_find_recoiler
    module function shower_generate_next_isr_branching &
         (shower) result (next_brancher)
      class(shower_t), intent(inout) :: shower
      type(parton_pointer_t) :: next_brancher
    end function shower_generate_next_isr_branching
    module subroutine shower_generate_fsr_for_partons_emitted_in_ISR (shower)
      class(shower_t), intent(inout) :: shower
    end subroutine shower_generate_fsr_for_partons_emitted_in_ISR
    module subroutine shower_execute_next_isr_branching (shower, prtp)
      class(shower_t), intent(inout) :: shower
      type(parton_pointer_t), intent(inout) :: prtp
    end subroutine shower_execute_next_isr_branching
    module function shower_get_ISR_scale (shower) result (scale)
      class(shower_t), intent(in) :: shower
      real(default) :: scale
    end function shower_get_ISR_scale
    module subroutine shower_set_max_isr_scale (shower, newscale)
      class(shower_t), intent(inout) :: shower
      real(default), intent(in) :: newscale
    end subroutine shower_set_max_isr_scale
    module subroutine shower_interaction_generate_fsr_2ton (shower, interaction)
      class(shower_t), intent(inout) :: shower
      type(shower_interaction_t), intent(inout) :: interaction
    end subroutine shower_interaction_generate_fsr_2ton
    recursive module subroutine shower_parton_pointer_array_generate_fsr_recursive &
         (shower, partons)
      class(shower_t), intent(inout) :: shower
      type(parton_pointer_t), dimension(:), allocatable, intent(inout) :: &
           partons
    end subroutine shower_parton_pointer_array_generate_fsr_recursive
    module subroutine shower_parton_pointer_array_generate_fsr &
         (shower, partons, partons_new)
      class(shower_t), intent(inout) :: shower
      type(parton_pointer_t), dimension(:), allocatable, intent(inout) :: &
           partons
      type(parton_pointer_t), dimension(:), allocatable, intent(out) :: &
           partons_new
    end subroutine shower_parton_pointer_array_generate_fsr
    module function shower_get_pdf (shower, mother, x, Q2, daughter) result (pdf)
      class(shower_t), intent(inout), target :: shower
      integer, intent(in) :: mother, daughter
      real(default), intent(in) :: x, Q2
      real(default) :: pdf
    end function shower_get_pdf
    module function shower_get_xpdf (shower, mother, x, Q2, daughter) result (pdf)
      class(shower_t), intent(inout), target :: shower
      integer, intent(in) :: mother, daughter
      real(default), intent(in) :: x, Q2
      real(default) :: pdf
    end function shower_get_xpdf
    module subroutine shower_converttopythia (shower)
      class(shower_t), intent(in) :: shower
    end subroutine shower_converttopythia
  end interface

end module shower_core
