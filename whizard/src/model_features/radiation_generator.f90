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

module radiation_generator

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use sorting, only: sort_abs
  use pdg_arrays
  use model_data
  use auto_components

  implicit none
  private

  public :: radiation_generator_t

  type :: pdg_sorter_t
     integer :: pdg
     logical :: checked = .false.
     integer :: associated_born = 0
  end type pdg_sorter_t

  type :: pdg_states_t
    type(pdg_array_t), dimension(:), allocatable :: pdg
    type(pdg_states_t), pointer :: next
    integer :: n_particles
  contains
    procedure :: init => pdg_states_init
    procedure :: add => pdg_states_add
    procedure :: get_n_states => pdg_states_get_n_states
  end type pdg_states_t

  type :: prt_queue_t
    type(string_t), dimension(:), allocatable :: prt_string
    type(prt_queue_t), pointer :: next => null ()
    type(prt_queue_t), pointer :: previous => null ()
    type(prt_queue_t), pointer :: front => null ()
    type(prt_queue_t), pointer :: current_prt => null ()
    type(prt_queue_t), pointer :: back => null ()
    integer :: n_lists = 0
  contains
    procedure :: null => prt_queue_null
    procedure :: append => prt_queue_append
    procedure :: get => prt_queue_get
    procedure :: get_last => prt_queue_get_last
    procedure :: reset => prt_queue_reset
    procedure :: check_for_same_prt_strings => prt_queue_check_for_same_prt_strings
    procedure :: contains => prt_queue_contains
    procedure :: write => prt_queue_write
  end type prt_queue_t

  type :: reshuffle_list_t
     integer, dimension(:), allocatable :: ii
     type(reshuffle_list_t), pointer :: next => null ()
  contains
    procedure :: write => reshuffle_list_write
    procedure :: append => reshuffle_list_append
    procedure :: is_empty => reshuffle_list_is_empty
    procedure :: get => reshuffle_list_get
    procedure :: reset => reshuffle_list_reset
  end type reshuffle_list_t

  type :: radiation_generator_t
    logical :: qcd_enabled = .false.
    logical :: qed_enabled = .false.
    logical :: is_gluon = .false.
    logical :: fs_gluon = .false.
    logical :: is_photon = .false.
    logical :: fs_photon = .false.
    logical :: only_final_state = .true.
    type(pdg_list_t) :: pl_in, pl_out
    type(pdg_list_t) :: pl_excluded_gauge_splittings
    type(split_constraints_t) :: constraints
    integer :: n_tot
    integer :: n_in, n_out
    integer :: n_loops
    integer :: n_light_quarks
    real(default) :: mass_sum
    type(prt_queue_t) :: prt_queue
    type(pdg_states_t) :: pdg_raw
    type(pdg_array_t), dimension(:), allocatable :: pdg_in_born, pdg_out_born
    type(if_table_t) :: if_table
    type(reshuffle_list_t) :: reshuffle_list
  contains
    generic :: init => init_pdg_list, init_pdg_array
    procedure :: init_pdg_list => radiation_generator_init_pdg_list
    procedure :: init_pdg_array => radiation_generator_init_pdg_array
    procedure :: set_initial_state_emissions => &
       radiation_generator_set_initial_state_emissions
    procedure :: setup_if_table => radiation_generator_setup_if_table
    generic :: reset_particle_content => reset_particle_content_pdg_array, &
                                         reset_particle_content_pdg_list
    procedure :: reset_particle_content_pdg_list => &
      radiation_generator_reset_particle_content_pdg_list
    procedure :: reset_particle_content_pdg_array => &
      radiation_generator_reset_particle_content_pdg_array
    procedure :: reset_reshuffle_list=> radiation_generator_reset_reshuffle_list
    procedure :: set_n => radiation_generator_set_n
    procedure :: set_constraints => radiation_generator_set_constraints
    procedure :: find_splittings => radiation_generator_find_splittings
    procedure :: generate_real_particle_strings &
         => radiation_generator_generate_real_particle_strings
    procedure :: contains_emissions => radiation_generator_contains_emissions
    procedure :: generate => radiation_generator_generate
    procedure :: generate_multiple => radiation_generator_generate_multiple
    procedure :: first_emission => radiation_generator_first_emission
    procedure :: append_emissions => radiation_generator_append_emissions
    procedure :: reset_queue => radiation_generator_reset_queue
    procedure :: get_n_gks_states => radiation_generator_get_n_gks_states
    procedure :: get_next_state => radiation_generator_get_next_state
    procedure :: get_emitter_indices => radiation_generator_get_emitter_indices
    procedure :: get_raw_states => radiation_generator_get_raw_states
    procedure :: save_born_raw => radiation_generator_save_born_raw
    procedure :: get_born_raw => radiation_generator_get_born_raw
  end type radiation_generator_t

  type :: prt_array_t
     type(string_t), dimension(:), allocatable :: prt
  end type prt_array_t

  type :: prt_table_t
     type(string_t), dimension(:), allocatable :: prt
  end type prt_table_t


  interface
    module subroutine pdg_states_init (states)
      class(pdg_states_t), intent(inout) :: states
    end subroutine pdg_states_init
    module subroutine pdg_states_add (states, pdg)
      class(pdg_states_t), intent(inout), target :: states
      type(pdg_array_t), dimension(:), intent(in) :: pdg
    end subroutine pdg_states_add
    module function pdg_states_get_n_states (states) result (n)
      class(pdg_states_t), intent(in), target :: states
      integer :: n
    end function pdg_states_get_n_states
    module subroutine prt_queue_null (queue)
      class(prt_queue_t), intent(out) :: queue
    end subroutine prt_queue_null
    module subroutine prt_queue_append (queue, prt_string)
      class(prt_queue_t), intent(inout) :: queue
      type(string_t), intent(in), dimension(:) :: prt_string
    end subroutine prt_queue_append
    module subroutine prt_queue_get (queue, prt_string)
      class(prt_queue_t), intent(inout) :: queue
      type(string_t), dimension(:), allocatable, intent(out) :: prt_string
    end subroutine prt_queue_get
    module subroutine prt_queue_get_last (queue, prt_string)
      class(prt_queue_t), intent(in) :: queue
      type(string_t), dimension(:), allocatable, intent(out) :: prt_string
    end subroutine prt_queue_get_last
    module subroutine prt_queue_reset (queue)
      class(prt_queue_t), intent(inout) :: queue
    end subroutine prt_queue_reset
    module function prt_queue_check_for_same_prt_strings (queue) result (val)
      class(prt_queue_t), intent(inout) :: queue
      logical :: val
    end function prt_queue_check_for_same_prt_strings
    module function prt_queue_contains (queue, prt_string) result (val)
      class(prt_queue_t), intent(in) :: queue
      type(string_t), intent(in), dimension(:) :: prt_string
      logical :: val
    end function prt_queue_contains
    module subroutine prt_queue_write (queue, unit)
      class(prt_queue_t), intent(in) :: queue
      integer, optional :: unit
    end subroutine prt_queue_write
    module subroutine reshuffle_list_write (rlist)
      class(reshuffle_list_t), intent(in) :: rlist
    end subroutine reshuffle_list_write
    module subroutine reshuffle_list_append (rlist, ii)
      class(reshuffle_list_t), intent(inout) :: rlist
      integer, dimension(:), allocatable, intent(in) :: ii
    end subroutine reshuffle_list_append
    elemental module function reshuffle_list_is_empty (rlist) result (is_empty)
      logical :: is_empty
      class(reshuffle_list_t), intent(in) :: rlist
    end function reshuffle_list_is_empty
    module function reshuffle_list_get (rlist, index) result (ii)
      integer, dimension(:), allocatable :: ii
      class(reshuffle_list_t), intent(inout) :: rlist
      integer, intent(in) :: index
    end function reshuffle_list_get
    module subroutine reshuffle_list_reset (rlist)
      class(reshuffle_list_t), intent(inout) :: rlist
    end subroutine reshuffle_list_reset
    module subroutine radiation_generator_init_pdg_list &
         (generator, pl_in, pl_out, pl_excluded_gauge_splittings, qcd, qed)
      class(radiation_generator_t), intent(inout) :: generator
      type(pdg_list_t), intent(in) :: pl_in, pl_out
      type(pdg_list_t), intent(in) :: pl_excluded_gauge_splittings
      logical, intent(in), optional :: qcd, qed
    end subroutine radiation_generator_init_pdg_list
    module subroutine radiation_generator_init_pdg_array &
         (generator, pdg_in, pdg_out, pdg_excluded_gauge_splittings, qcd, qed)
      class(radiation_generator_t), intent(inout) :: generator
      type(pdg_array_t), intent(in), dimension(:) :: pdg_in, pdg_out
      type(pdg_array_t), intent(in), dimension(:) :: pdg_excluded_gauge_splittings
      logical, intent(in), optional :: qcd, qed
    end subroutine radiation_generator_init_pdg_array
    module subroutine radiation_generator_set_initial_state_emissions (generator)
       class(radiation_generator_t), intent(inout) :: generator
    end subroutine radiation_generator_set_initial_state_emissions
    module subroutine radiation_generator_setup_if_table (generator, model)
      class(radiation_generator_t), intent(inout) :: generator
      class(model_data_t), intent(in), target :: model
    end subroutine radiation_generator_setup_if_table
    module subroutine radiation_generator_reset_particle_content_pdg_list (generator, pl)
      class(radiation_generator_t), intent(inout) :: generator
      type(pdg_list_t), intent(in) :: pl
    end subroutine radiation_generator_reset_particle_content_pdg_list
    module subroutine radiation_generator_reset_particle_content_pdg_array (generator, pdg)
      class(radiation_generator_t), intent(inout) :: generator
      type(pdg_array_t), intent(in), dimension(:) :: pdg
    end subroutine radiation_generator_reset_particle_content_pdg_array
    module subroutine radiation_generator_reset_reshuffle_list (generator)
      class(radiation_generator_t), intent(inout) :: generator
    end subroutine radiation_generator_reset_reshuffle_list
    module subroutine radiation_generator_set_n (generator, n_in, n_out, n_loops)
      class(radiation_generator_t), intent(inout) :: generator
      integer, intent(in) :: n_in, n_out, n_loops
    end subroutine radiation_generator_set_n
    module subroutine radiation_generator_set_constraints &
         (generator, set_n_loop, set_mass_sum, &
          set_selected_particles, set_required_particles)
      class(radiation_generator_t), intent(inout), target :: generator
      logical, intent(in) :: set_n_loop
      logical, intent(in) :: set_mass_sum
      logical, intent(in) :: set_selected_particles
      logical, intent(in) :: set_required_particles
    end subroutine radiation_generator_set_constraints
    module subroutine radiation_generator_find_splittings (generator)
      class(radiation_generator_t), intent(inout) :: generator
    end subroutine radiation_generator_find_splittings
    module subroutine radiation_generator_generate_real_particle_strings &
         (generator, prt_tot_in, prt_tot_out)
      class(radiation_generator_t), intent(inout) :: generator
      type(string_t), intent(out), dimension(:), allocatable :: prt_tot_in, prt_tot_out
    end subroutine radiation_generator_generate_real_particle_strings
    module function radiation_generator_contains_emissions (generator) result (has_em)
      logical :: has_em
      class(radiation_generator_t), intent(in) :: generator
    end function radiation_generator_contains_emissions
    module subroutine radiation_generator_generate (generator, prt_in, prt_out)
      class(radiation_generator_t), intent(inout) :: generator
      type(string_t), intent(out), dimension(:), allocatable :: prt_in, prt_out
    end subroutine radiation_generator_generate
    module subroutine radiation_generator_generate_multiple &
         (generator, max_multiplicity, model)
      class(radiation_generator_t), intent(inout) :: generator
      integer, intent(in) :: max_multiplicity
      class(model_data_t), intent(in), target :: model
    end subroutine radiation_generator_generate_multiple
    module subroutine radiation_generator_first_emission (generator, model)
      class(radiation_generator_t), intent(inout) :: generator
      class(model_data_t), intent(in), target :: model
    end subroutine radiation_generator_first_emission
    module subroutine radiation_generator_append_emissions &
         (generator, max_multiplicity, model)
      class(radiation_generator_t), intent(inout) :: generator
      integer, intent(in) :: max_multiplicity
      class(model_data_t), intent(in), target :: model
    end subroutine radiation_generator_append_emissions
    module subroutine radiation_generator_reset_queue (generator)
      class(radiation_generator_t), intent(inout) :: generator
    end subroutine radiation_generator_reset_queue
    module function radiation_generator_get_n_gks_states (generator) result (n)
      class(radiation_generator_t), intent(in) :: generator
      integer :: n
    end function radiation_generator_get_n_gks_states
    module function radiation_generator_get_next_state (generator) result (prt_string)
      class(radiation_generator_t), intent(inout) :: generator
      type(string_t), dimension(:), allocatable :: prt_string
    end function radiation_generator_get_next_state
    module subroutine radiation_generator_get_emitter_indices (generator, indices)
      class(radiation_generator_t), intent(in) :: generator
      integer, dimension(:), allocatable, intent(out) :: indices
    end subroutine radiation_generator_get_emitter_indices
    module function radiation_generator_get_raw_states (generator) result (raw_states)
      class(radiation_generator_t), intent(in), target :: generator
      integer, dimension(:,:), allocatable :: raw_states
    end function radiation_generator_get_raw_states
    module subroutine radiation_generator_save_born_raw (generator, pdg_in, pdg_out)
      class(radiation_generator_t), intent(inout) :: generator
      type(pdg_array_t), dimension(:), allocatable, intent(in) :: pdg_in, pdg_out
    end subroutine radiation_generator_save_born_raw
    module function radiation_generator_get_born_raw (generator) result (flv_born)
      class(radiation_generator_t), intent(in) :: generator
      integer, dimension(:,:), allocatable :: flv_born
    end function radiation_generator_get_born_raw
  end interface

end module radiation_generator
