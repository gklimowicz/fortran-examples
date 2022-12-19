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

module cascades

  use kinds, only: default
  use kinds, only: TC, i8, i32
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use physics_defs, only: SCALAR, SPINOR, VECTOR, VECTORSPINOR, TENSOR
  use physics_defs, only: UNDEFINED
  use model_data
  use flavors

  use resonances, only: resonance_info_t
  use resonances, only: resonance_history_t
  use resonances, only: resonance_history_set_t
  use phs_forests

  implicit none
  private

  public :: hash_entry_init
  public :: cascade_set_t
  interface cascade_set_init
     module procedure cascade_set_init_base
     module procedure cascade_set_init_from_cascade
  end interface
  public :: cascade_set_is_valid
  public :: cascade_set_final
  public :: cascade_set_write_process_bincode_format
  public :: cascade_set_write_file_format
  public :: cascade_set_write_graph_format
  public :: cascade_set_write
  public :: cascade_set_generate
  public :: phase_space_vanishes
  public :: assignment(=)
  public :: cascade_set_get_n_trees
  public :: cascade_set_get_resonance_histories

  integer, parameter :: &
       & EXTERNAL_PRT = -1, &
       & NO_MAPPING = 0, S_CHANNEL = 1, T_CHANNEL =  2, U_CHANNEL = 3, &
       & RADIATION = 4, COLLINEAR = 5, INFRARED = 6, &
       & STEP_MAPPING_E = 11, STEP_MAPPING_H = 12, &
       & ON_SHELL = 99
  real, parameter, public :: CASCADE_SET_FILL_RATIO = 0.1
  integer, parameter, public :: MAX_WARN_RESONANCE = 50

  type :: cascade_t
     private
     ! counters
     integer :: index = 0
     integer :: grove = 0
     ! status
     logical :: active = .false.
     logical :: complete = .false.
     logical :: incoming = .false.
     ! this node
     integer(TC) :: bincode = 0
     type(flavor_t) :: flv
     integer :: pdg = UNDEFINED
     logical :: is_vector = .false.
     real(default) :: m_min = 0
     real(default) :: m_rea = 0
     real(default) :: m_eff = 0
     integer :: mapping = NO_MAPPING
     logical :: on_shell = .false.
     logical :: resonant = .false.
     logical :: log_enhanced = .false.
     logical :: t_channel = .false.
     ! global tree properties
     integer :: multiplicity = 0
     integer :: internal = 0
     integer :: n_off_shell = 0
     integer :: n_resonances = 0
     integer :: n_log_enhanced = 0
     integer :: n_t_channel = 0
     integer :: res_hash = 0
     ! the sub-node tree
     integer :: depth = 0
     integer(TC), dimension(:), allocatable :: tree
     integer, dimension(:), allocatable :: tree_pdg
     integer, dimension(:), allocatable :: tree_mapping
     logical, dimension(:), allocatable :: tree_resonant
     ! branch connections
     logical :: has_children = .false.
     type(cascade_t), pointer :: daughter1 => null ()
     type(cascade_t), pointer :: daughter2 => null ()
     type(cascade_t), pointer :: mother => null ()
     ! next in list
     type(cascade_t), pointer :: next => null ()
   contains
     procedure :: extract_resonance_history => cascade_extract_resonance_history
  end type cascade_t

  type :: cascade_p
     type(cascade_t), pointer :: cascade => null ()
     type(cascade_p), pointer :: next => null ()
  end type cascade_p

  type :: hash_entry_t
     integer(i32) :: hashval = 0
     integer(i8), dimension(:), allocatable :: key
     type(cascade_p), pointer :: first => null ()
     type(cascade_p), pointer :: last => null ()
  end type hash_entry_t

  type :: cascade_set_t
     private
     class(model_data_t), pointer :: model
     integer :: n_in, n_out, n_tot
     type(flavor_t), dimension(:,:), allocatable :: flv
     integer :: depth_out, depth_tot
     real(default) :: sqrts = 0
     real(default) :: m_threshold_s = 0
     real(default) :: m_threshold_t = 0
     integer :: off_shell = 0
     integer :: t_channel = 0
     logical :: keep_nonresonant
     integer :: n_groves = 0
     ! The cascade list
     type(cascade_t), pointer :: first => null ()
     type(cascade_t), pointer :: last => null ()
     type(cascade_t), pointer :: first_t => null ()
     type(cascade_t), pointer :: first_k => null ()
     ! The hashtable
     integer :: n_entries = 0
     real :: fill_ratio = 0
     integer :: n_entries_max = 0
     integer(i32) :: mask = 0
     logical :: fatal_beam_decay = .true.
     type(hash_entry_t), dimension(:), allocatable :: entry
  end type cascade_set_t


  interface operator(.disjunct.)
     module procedure cascade_disjunct
  end interface

  interface operator(.match.)
     module procedure pdg_match
  end interface
  interface cascade_set_add_outgoing
     module procedure cascade_set_add_outgoing1
     module procedure cascade_set_add_outgoing2
  end interface

  interface cascade_set_add_incoming
     module procedure cascade_set_add_incoming0
     module procedure cascade_set_add_incoming1
  end interface


  interface
    module function cascade_disjunct (cascade1, cascade2) result (flag)
      logical :: flag
      type(cascade_t), intent(in) :: cascade1, cascade2
    end function cascade_disjunct
    module subroutine hash_entry_init (entry, entry_in)
      type(hash_entry_t), intent(out) :: entry
      type(hash_entry_t), intent(in) :: entry_in
    end subroutine hash_entry_init
    elemental module function pdg_match (pdg1, pdg2) result (flag)
      logical :: flag
      integer(TC), intent(in) :: pdg1, pdg2
    end function pdg_match
    module subroutine cascade_set_init_from_cascade &
         (cascade_set, cascade_set_in)
      type(cascade_set_t), intent(out) :: cascade_set
      type(cascade_set_t), intent(in), target :: cascade_set_in
    end subroutine cascade_set_init_from_cascade
    module function cascade_set_is_valid (cascade_set) result (flag)
      logical :: flag
      type(cascade_set_t), intent(in) :: cascade_set
    end function cascade_set_is_valid
    module subroutine cascade_set_init_base (cascade_set, model, &
         n_in, n_out, phs_par, fatal_beam_decay, flv)
      type(cascade_set_t), intent(out) :: cascade_set
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_in, n_out
      type(phs_parameters_t), intent(in) :: phs_par
      logical, intent(in) :: fatal_beam_decay
      type(flavor_t), dimension(:,:), intent(in), optional :: flv
    end subroutine cascade_set_init_base
    module subroutine cascade_set_final (cascade_set)
      type(cascade_set_t), intent(inout), target :: cascade_set
    end subroutine cascade_set_final
    module subroutine cascade_set_write_process_bincode_format &
         (cascade_set, unit)
      type(cascade_set_t), intent(in), target :: cascade_set
      integer, intent(in), optional :: unit
    end subroutine cascade_set_write_process_bincode_format
    module subroutine cascade_set_write_file_format (cascade_set, unit)
      type(cascade_set_t), intent(in), target :: cascade_set
      integer, intent(in), optional :: unit
    end subroutine cascade_set_write_file_format
    module subroutine cascade_set_write_graph_format &
        (cascade_set, filename, process_id, unit)
      type(cascade_set_t), intent(in), target :: cascade_set
      type(string_t), intent(in) :: filename, process_id
      integer, intent(in), optional :: unit
    end subroutine cascade_set_write_graph_format
    module subroutine cascade_set_write &
         (cascade_set, unit, active_only, complete_only)
      type(cascade_set_t), intent(in), target :: cascade_set
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: active_only, complete_only
    end subroutine cascade_set_write
    module subroutine cascade_set_add_outgoing1 (cascade_set, flv)
      type(cascade_set_t), intent(inout), target :: cascade_set
      type(flavor_t), dimension(:), intent(in) :: flv
    end subroutine cascade_set_add_outgoing1
    module subroutine cascade_set_add_outgoing2 (cascade_set, flv)
      type(cascade_set_t), intent(inout), target :: cascade_set
      type(flavor_t), dimension(:,:), intent(in) :: flv
    end subroutine cascade_set_add_outgoing2
    module subroutine cascade_set_add_incoming1 (cascade_set, n1, n2, pos, flv)
      type(cascade_set_t), intent(inout), target :: cascade_set
      integer, intent(out) :: n1, n2
      integer, intent(in) :: pos
      type(flavor_t), dimension(:), intent(in) :: flv
    end subroutine cascade_set_add_incoming1
    module subroutine cascade_set_add_incoming0 (cascade_set, n1, n2, pos, flv)
      type(cascade_set_t), intent(inout), target :: cascade_set
      integer, intent(out) :: n1, n2
      integer, intent(in) :: pos
      type(flavor_t), intent(in) :: flv
    end subroutine cascade_set_add_incoming0
    module subroutine cascade_set_generate &
         (cascade_set, model, n_in, n_out, flv, phs_par, fatal_beam_decay)
      type(cascade_set_t), intent(out) :: cascade_set
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_in, n_out
      type(flavor_t), dimension(:,:), intent(in) :: flv
      type(phs_parameters_t), intent(in) :: phs_par
      logical, intent(in) :: fatal_beam_decay
    end subroutine cascade_set_generate
    module function phase_space_vanishes (sqrts, n_in, flv) result (flag)
      logical :: flag
      real(default), intent(in) :: sqrts
      integer, intent(in) :: n_in
      type(flavor_t), dimension(:,:), intent(in) :: flv
    end function phase_space_vanishes
    module subroutine cascade_extract_resonance_history &
         (cascade, res_hist, model, n_out)
      class(cascade_t), intent(in), target :: cascade
      type(resonance_history_t), intent(out) :: res_hist
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: n_out
    end subroutine cascade_extract_resonance_history
    module function cascade_set_get_n_trees (cascade_set) result (n)
      type(cascade_set_t), intent(in), target :: cascade_set
      integer :: n
    end function cascade_set_get_n_trees
    module subroutine cascade_set_get_resonance_histories &
         (cascade_set, n_filter, res_hists)
      type(cascade_set_t), intent(in), target :: cascade_set
      integer, intent(in), optional :: n_filter
      type(resonance_history_t), dimension(:), allocatable, intent(out) :: &
           res_hists
    end subroutine cascade_set_get_resonance_histories
  end interface

end module cascades
