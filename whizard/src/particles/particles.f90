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

module particles

  use kinds, only: default, double
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use lorentz
  use phs_points, only: phs_point_t, assignment(=)
  use model_data
  use flavors
  use colors
  use helicities
  use quantum_numbers
  use state_matrices
  use interactions
  use subevents
  use polarizations

  implicit none
  private

  public :: particle_t
  public :: particle_set_t
  public :: pacify

  integer, parameter, public :: PRT_UNPOLARIZED = 0
  integer, parameter, public :: PRT_DEFINITE_HELICITY = 1
  integer, parameter, public :: PRT_GENERIC_POLARIZATION = 2


  type :: particle_t
     !private
     integer :: status = PRT_UNDEFINED
     integer :: polarization = PRT_UNPOLARIZED
     type(flavor_t) :: flv
     type(color_t) :: col
     type(helicity_t) :: hel
     type(polarization_t) :: pol
     type(vector4_t) :: p = vector4_null
     real(default) :: p2 = 0
     type(vector4_t), allocatable :: vertex
     real(default), allocatable :: lifetime
     integer, dimension(:), allocatable :: parent
     integer, dimension(:), allocatable :: child
   contains
    generic :: init => init_particle
    procedure :: init_particle => particle_init_particle
    generic :: init => init_external
    procedure :: init_external => particle_init_external
    generic :: init => init_state
    procedure :: init_state => particle_init_state
    procedure :: final => particle_final
    procedure :: write => particle_write
    procedure :: write_raw => particle_write_raw
    procedure :: read_raw => particle_read_raw
    procedure :: reset_status => particle_reset_status
    procedure :: set_color => particle_set_color
    procedure :: set_flavor => particle_set_flavor
    procedure :: set_helicity => particle_set_helicity
    procedure :: set_pol => particle_set_pol
    procedure :: set_model => particle_set_model
    procedure :: set_momentum => particle_set_momentum
    procedure :: set_resonance_flag => particle_set_resonance_flag
    procedure :: set_children => particle_set_children
    procedure :: set_parents => particle_set_parents
    procedure :: add_child => particle_add_child
    procedure :: add_children => particle_add_children
    procedure :: set_status => particle_set_status
    procedure :: set_polarization => particle_set_polarization
    generic :: set_vertex => set_vertex_from_vector3, set_vertex_from_xyz, &
         set_vertex_from_vector4, set_vertex_from_xyzt
    procedure :: set_vertex_from_vector4 => particle_set_vertex_from_vector4
    procedure :: set_vertex_from_vector3 => particle_set_vertex_from_vector3
    procedure :: set_vertex_from_xyzt => particle_set_vertex_from_xyzt
    procedure :: set_vertex_from_xyz => particle_set_vertex_from_xyz
    procedure :: set_lifetime => particle_set_lifetime
    procedure :: get_status => particle_get_status
    procedure :: is_real => particle_is_real
    procedure :: is_colored => particle_is_colored
    procedure :: is_hadronic_beam_remnant => particle_is_hadronic_beam_remnant
    procedure :: is_beam_remnant => particle_is_beam_remnant
    procedure :: get_polarization_status => particle_get_polarization_status
    procedure :: get_pdg => particle_get_pdg
    procedure :: get_color => particle_get_color
    procedure :: get_polarization => particle_get_polarization
    procedure :: get_flv => particle_get_flv
    procedure :: get_col => particle_get_col
    procedure :: get_hel => particle_get_hel
    procedure :: get_helicity => particle_get_helicity
    procedure :: get_n_parents => particle_get_n_parents
    procedure :: get_n_children => particle_get_n_children
    procedure :: get_parents => particle_get_parents
    procedure :: get_children => particle_get_children
    procedure :: has_children => particle_has_children
    procedure :: has_parents => particle_has_parents
    procedure :: get_momentum => particle_get_momentum
    procedure :: get_p2 => particle_get_p2
    procedure :: get_vertex => particle_get_vertex
    procedure :: get_lifetime => particle_get_lifetime
    procedure :: momentum_to_pythia6 => particle_momentum_to_pythia6
  end type particle_t

  type :: particle_set_t
     ! private !!!
     integer :: n_beam = 0
     integer :: n_in  = 0
     integer :: n_vir = 0
     integer :: n_out = 0
     integer :: n_tot = 0
     integer :: factorization_mode = FM_IGNORE_HELICITY
     type(particle_t), dimension(:), allocatable :: prt
     type(state_matrix_t) :: correlated_state
   contains
     generic :: init => init_interaction
     procedure :: init_interaction => particle_set_init_interaction
     generic :: assignment(=) => init_particle_set
     generic :: init => init_particle_set
     procedure :: init_particle_set => particle_set_init_particle_set
     procedure :: set_model => particle_set_set_model
     procedure :: final => particle_set_final
     procedure :: basic_init => particle_set_basic_init
     procedure :: init_direct => particle_set_init_direct
     procedure :: transfer => particle_set_transfer
     procedure :: insert => particle_set_insert
     procedure :: recover_color => particle_set_recover_color
     generic :: get_color => get_color_all
     generic :: get_color => get_color_indices
     procedure :: get_color_all => particle_set_get_color_all
     procedure :: get_color_indices => particle_set_get_color_indices
     generic :: set_color => set_color_single
     generic :: set_color => set_color_indices
     generic :: set_color => set_color_all
     procedure :: set_color_single => particle_set_set_color_single
     procedure :: set_color_indices => particle_set_set_color_indices
     procedure :: set_color_all => particle_set_set_color_all
     procedure :: find_prt_invalid_color => particle_set_find_prt_invalid_color
     generic :: get_momenta => get_momenta_all
     generic :: get_momenta => get_momenta_indices
     procedure :: get_momenta_all => particle_set_get_momenta_all
     procedure :: get_momenta_indices => particle_set_get_momenta_indices
     generic :: set_momentum => set_momentum_single
     generic :: set_momentum => set_momentum_indices
     generic :: set_momentum => set_momentum_all
     procedure :: set_momentum_single => particle_set_set_momentum_single
     procedure :: set_momentum_indices => particle_set_set_momentum_indices
     procedure :: set_momentum_all => particle_set_set_momentum_all
     procedure :: recover_momentum => particle_set_recover_momentum
     procedure :: replace_incoming_momenta => particle_set_replace_incoming_momenta
     procedure :: replace_outgoing_momenta => particle_set_replace_outgoing_momenta
     procedure :: get_outgoing_momenta => particle_set_get_outgoing_momenta
     procedure :: parent_add_child => particle_set_parent_add_child
     procedure :: build_radiation => particle_set_build_radiation
     procedure :: write => particle_set_write
     procedure :: write_raw => particle_set_write_raw
     procedure :: read_raw => particle_set_read_raw
     procedure :: get_real_parents => particle_set_get_real_parents
     procedure :: get_real_children => particle_set_get_real_children
     procedure :: get_n_beam => particle_set_get_n_beam
     procedure :: get_n_in => particle_set_get_n_in
     procedure :: get_n_vir => particle_set_get_n_vir
     procedure :: get_n_out => particle_set_get_n_out
     procedure :: get_n_tot => particle_set_get_n_tot
     procedure :: get_n_remnants => particle_set_get_n_remnants
     procedure :: get_particle => particle_set_get_particle
     procedure :: get_indices => particle_set_get_indices
     procedure :: get_in_and_out_momenta => particle_set_get_in_and_out_momenta
     procedure :: without_hadronic_remnants => &
          particle_set_without_hadronic_remnants
     procedure :: without_remnants => particle_set_without_remnants
     procedure :: find_particle => particle_set_find_particle
     procedure :: reverse_find_particle => particle_set_reverse_find_particle
     procedure :: remove_duplicates => particle_set_remove_duplicates
     procedure :: reset_status => particle_set_reset_status
     procedure :: reduce => particle_set_reduce
     procedure :: filter_particles => particle_set_filter_particles
     procedure :: to_hepevt_form => particle_set_to_hepevt_form
     procedure :: fill_interaction => particle_set_fill_interaction
     procedure :: assign_vertices => particle_set_assign_vertices
     procedure :: to_subevt => particle_set_to_subevt
     procedure :: replace => particle_set_replace
     procedure :: order_color_lines => particle_set_order_color_lines
  end type particle_set_t

  type :: particle_entry_t
     integer :: src = 0
     integer :: status = 0
     integer :: orig = 0
     integer :: copy = 0
  end type particle_entry_t


  interface pacify
     module procedure pacify_particle
     module procedure pacify_particle_set
  end interface pacify


  interface
    module subroutine particle_init_particle (prt_out, prt_in)
      class(particle_t), intent(out) :: prt_out
      type(particle_t), intent(in) :: prt_in
    end subroutine particle_init_particle
    module subroutine particle_init_external &
           (particle, status, pdg, model, col, anti_col, mom)
      class(particle_t), intent(out) :: particle
      integer, intent(in) :: status, pdg, col, anti_col
      class(model_data_t), pointer, intent(in) :: model
      type(vector4_t), intent(in) :: mom
    end subroutine particle_init_external
    module subroutine particle_init_state (prt, state, status, mode)
      class(particle_t), intent(out) :: prt
      type(state_matrix_t), intent(in), target :: state
      integer, intent(in) :: status, mode
    end subroutine particle_init_state
    module subroutine particle_final (prt)
      class(particle_t), intent(inout) :: prt
    end subroutine particle_final
    module subroutine particle_write (prt, unit, testflag, compressed, polarization)
      class(particle_t), intent(in) :: prt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag, compressed, polarization
    end subroutine particle_write
    module subroutine particle_write_raw (prt, u)
      class(particle_t), intent(in) :: prt
      integer, intent(in) :: u
    end subroutine particle_write_raw
    module subroutine particle_read_raw (prt, u, iostat)
      class(particle_t), intent(out) :: prt
      integer, intent(in) :: u
      integer, intent(out) :: iostat
    end subroutine particle_read_raw
    elemental module subroutine particle_reset_status (prt, status)
      class(particle_t), intent(inout) :: prt
      integer, intent(in) :: status
    end subroutine particle_reset_status
    elemental module subroutine particle_set_color (prt, col)
      class(particle_t), intent(inout) :: prt
      type(color_t), intent(in) :: col
    end subroutine particle_set_color
    module subroutine particle_set_flavor (prt, flv)
      class(particle_t), intent(inout) :: prt
      type(flavor_t), intent(in) :: flv
    end subroutine particle_set_flavor
    module subroutine particle_set_helicity (prt, hel)
      class(particle_t), intent(inout) :: prt
      type(helicity_t), intent(in) :: hel
    end subroutine particle_set_helicity
    module subroutine particle_set_pol (prt, pol)
      class(particle_t), intent(inout) :: prt
      type(polarization_t), intent(in) :: pol
    end subroutine particle_set_pol
    module subroutine particle_set_model (prt, model)
      class(particle_t), intent(inout) :: prt
      class(model_data_t), intent(in), target :: model
    end subroutine particle_set_model
    elemental module subroutine particle_set_momentum (prt, p, p2, on_shell)
      class(particle_t), intent(inout) :: prt
      type(vector4_t), intent(in) :: p
      real(default), intent(in), optional :: p2
      logical, intent(in), optional :: on_shell
    end subroutine particle_set_momentum
    elemental module subroutine particle_set_resonance_flag (prt, resonant)
      class(particle_t), intent(inout) :: prt
      logical, intent(in) :: resonant
    end subroutine particle_set_resonance_flag
    module subroutine particle_set_children (prt, idx)
      class(particle_t), intent(inout) :: prt
      integer, dimension(:), intent(in) :: idx
    end subroutine particle_set_children
    module subroutine particle_set_parents (prt, idx)
      class(particle_t), intent(inout) :: prt
      integer, dimension(:), intent(in) :: idx
    end subroutine particle_set_parents
    module subroutine particle_add_child (prt, new_child)
      class(particle_t), intent(inout) :: prt
      integer, intent(in) :: new_child
    end subroutine particle_add_child
    module subroutine particle_add_children (prt, new_child)
      class(particle_t), intent(inout) :: prt
      integer, dimension(:), intent(in) :: new_child
    end subroutine particle_add_children
    elemental module subroutine particle_set_status (prt, status)
      class(particle_t), intent(inout) :: prt
      integer, intent(in) :: status
    end subroutine particle_set_status
    module subroutine particle_set_polarization (prt, polarization)
      class(particle_t), intent(inout) :: prt
      integer, intent(in) :: polarization
    end subroutine particle_set_polarization
    module subroutine particle_set_vertex_from_vector4 (prt, vertex)
      class(particle_t), intent(inout) :: prt
      type(vector4_t), intent(in) :: vertex
    end subroutine particle_set_vertex_from_vector4
    module subroutine particle_set_vertex_from_vector3 (prt, vertex)
      class(particle_t), intent(inout) :: prt
      type(vector3_t), intent(in) :: vertex
    end subroutine particle_set_vertex_from_vector3
    module subroutine particle_set_vertex_from_xyzt (prt, vx, vy, vz, t)
      class(particle_t), intent(inout) :: prt
      real(default), intent(in) :: vx, vy, vz, t
    end subroutine particle_set_vertex_from_xyzt
    module subroutine particle_set_vertex_from_xyz (prt, vx, vy, vz)
      class(particle_t), intent(inout) :: prt
      real(default), intent(in) :: vx, vy, vz
    end subroutine particle_set_vertex_from_xyz
    elemental module subroutine particle_set_lifetime (prt, lifetime)
      class(particle_t), intent(inout) :: prt
      real(default), intent(in) :: lifetime
    end subroutine particle_set_lifetime
    elemental module function particle_get_status (prt) result (status)
      integer :: status
      class(particle_t), intent(in) :: prt
    end function particle_get_status
    elemental module function particle_is_real (prt, keep_beams) result (flag)
      logical :: flag, kb
      class(particle_t), intent(in) :: prt
      logical, intent(in), optional :: keep_beams
    end function particle_is_real
    elemental module function particle_is_colored (particle) result (flag)
      logical :: flag
      class(particle_t), intent(in) :: particle
    end function particle_is_colored
    elemental module function particle_is_hadronic_beam_remnant (particle) result (flag)
      class(particle_t), intent(in) :: particle
      logical :: flag
    end function particle_is_hadronic_beam_remnant
    elemental module function particle_is_beam_remnant (particle) result (flag)
      class(particle_t), intent(in) :: particle
      logical :: flag
    end function particle_is_beam_remnant
    elemental module function particle_get_polarization_status (prt) result (status)
      integer :: status
      class(particle_t), intent(in) :: prt
    end function particle_get_polarization_status
    elemental module function particle_get_pdg (prt) result (pdg)
      integer :: pdg
      class(particle_t), intent(in) :: prt
    end function particle_get_pdg
    pure module function particle_get_color (prt) result (col)
      integer, dimension(2) :: col
      class(particle_t), intent(in) :: prt
    end function particle_get_color
    module function particle_get_polarization (prt) result (pol)
      class(particle_t), intent(in) :: prt
      type(polarization_t) :: pol
    end function particle_get_polarization
    module function particle_get_flv (prt) result (flv)
      class(particle_t), intent(in) :: prt
      type(flavor_t) :: flv
    end function particle_get_flv
    module function particle_get_col (prt) result (col)
      class(particle_t), intent(in) :: prt
      type(color_t) :: col
    end function particle_get_col
    module function particle_get_hel (prt) result (hel)
      class(particle_t), intent(in) :: prt
      type(helicity_t) :: hel
    end function particle_get_hel
    elemental module function particle_get_helicity (prt) result (hel)
      integer :: hel
      class(particle_t), intent(in) :: prt
    end function particle_get_helicity
    elemental module function particle_get_n_parents (prt) result (n)
      integer :: n
      class(particle_t), intent(in) :: prt
    end function particle_get_n_parents
    elemental module function particle_get_n_children (prt) result (n)
      integer :: n
      class(particle_t), intent(in) :: prt
    end function particle_get_n_children
    module function particle_get_parents (prt) result (parent)
      class(particle_t), intent(in) :: prt
      integer, dimension(:), allocatable :: parent
    end function particle_get_parents
    module function particle_get_children (prt) result (child)
      class(particle_t), intent(in) :: prt
      integer, dimension(:), allocatable :: child
    end function particle_get_children
    elemental module function particle_has_children (prt) result (has_children)
      logical :: has_children
      class(particle_t), intent(in) :: prt
    end function particle_has_children
    elemental module function particle_has_parents (prt) result (has_parents)
      logical :: has_parents
      class(particle_t), intent(in) :: prt
    end function particle_has_parents
    elemental module function particle_get_momentum (prt) result (p)
      type(vector4_t) :: p
      class(particle_t), intent(in) :: prt
    end function particle_get_momentum
    elemental module function particle_get_p2 (prt) result (p2)
      real(default) :: p2
      class(particle_t), intent(in) :: prt
    end function particle_get_p2
    elemental module function particle_get_vertex (prt) result (vtx)
      type(vector4_t) :: vtx
      class(particle_t), intent(in) :: prt
    end function particle_get_vertex
    elemental module function particle_get_lifetime (prt) result (lifetime)
      real(default) :: lifetime
      class(particle_t), intent(in) :: prt
    end function particle_get_lifetime
    pure module function particle_momentum_to_pythia6 (prt) result (p)
      real(double), dimension(1:5) :: p
      class(particle_t), intent(in) :: prt
    end function particle_momentum_to_pythia6
    module subroutine particle_set_init_interaction &
         (particle_set, is_valid, int, int_flows, mode, x, &
          keep_correlations, keep_virtual, n_incoming, qn_select)
      class(particle_set_t), intent(out) :: particle_set
      logical, intent(out) :: is_valid
      type(interaction_t), intent(in), target :: int, int_flows
      integer, intent(in) :: mode
      real(default), dimension(2), intent(in) :: x
      logical, intent(in) :: keep_correlations, keep_virtual
      integer, intent(in), optional :: n_incoming
      type(quantum_numbers_t), dimension(:), intent(in), optional :: qn_select
    end subroutine particle_set_init_interaction
    module subroutine particle_set_init_particle_set (pset_out, pset_in)
      class(particle_set_t), intent(out) :: pset_out
      type(particle_set_t), intent(in) :: pset_in
    end subroutine particle_set_init_particle_set
    module subroutine particle_set_set_model (particle_set, model)
      class(particle_set_t), intent(inout) :: particle_set
      class(model_data_t), intent(in), target :: model
    end subroutine particle_set_set_model
    module subroutine particle_set_final (particle_set)
      class(particle_set_t), intent(inout) :: particle_set
    end subroutine particle_set_final
    module subroutine particle_set_basic_init (particle_set, n_beam, n_in, n_vir, n_out)
      class(particle_set_t), intent(out) :: particle_set
      integer, intent(in) :: n_beam, n_in, n_vir, n_out
    end subroutine particle_set_basic_init
    module subroutine particle_set_init_direct (particle_set, &
         n_beam, n_in, n_rem, n_vir, n_out, pdg, model)
      class(particle_set_t), intent(out) :: particle_set
      integer, intent(in) :: n_beam
      integer, intent(in) :: n_in
      integer, intent(in) :: n_rem
      integer, intent(in) :: n_vir
      integer, intent(in) :: n_out
      integer, dimension(:), intent(in) :: pdg
      class(model_data_t), intent(in), target :: model
    end subroutine particle_set_init_direct
    module subroutine particle_set_transfer (pset, source, n_new, map)
      class(particle_set_t), intent(out) :: pset
      class(particle_set_t), intent(in) :: source
      integer, intent(in) :: n_new
      integer, dimension(:), intent(in) :: map
    end subroutine particle_set_transfer
    module subroutine particle_set_insert (pset, i, status, flv, child)
      class(particle_set_t), intent(inout) :: pset
      integer, intent(in) :: i
      integer, intent(in) :: status
      type(flavor_t), intent(in) :: flv
      integer, dimension(:), intent(in) :: child
    end subroutine particle_set_insert
    module subroutine particle_set_recover_color (pset, i)
      class(particle_set_t), intent(inout) :: pset
      integer, intent(in) :: i
    end subroutine particle_set_recover_color
    module function particle_set_get_color_all (particle_set) result (col)
      class(particle_set_t), intent(in) :: particle_set
      type(color_t), dimension(:), allocatable :: col
    end function particle_set_get_color_all
    module function particle_set_get_color_indices (particle_set, indices) result (col)
      type(color_t), dimension(:), allocatable :: col
      class(particle_set_t), intent(in) :: particle_set
      integer, intent(in), dimension(:), allocatable :: indices
    end function particle_set_get_color_indices
    module subroutine particle_set_set_color_single (particle_set, i, col)
      class(particle_set_t), intent(inout) :: particle_set
      integer, intent(in) :: i
      type(color_t), intent(in) :: col
    end subroutine particle_set_set_color_single
    module subroutine particle_set_set_color_indices (particle_set, indices, col)
      class(particle_set_t), intent(inout) :: particle_set
      integer, dimension(:), intent(in) :: indices
      type(color_t), dimension(:), intent(in) :: col
    end subroutine particle_set_set_color_indices
    module subroutine particle_set_set_color_all (particle_set, col)
      class(particle_set_t), intent(inout) :: particle_set
      type(color_t), dimension(:), intent(in) :: col
    end subroutine particle_set_set_color_all
    module subroutine particle_set_find_prt_invalid_color (particle_set, index, prt)
      class(particle_set_t), intent(in) :: particle_set
      integer, dimension(:), allocatable, intent(out) :: index
      type(particle_t), dimension(:), allocatable, intent(out), optional :: prt
    end subroutine particle_set_find_prt_invalid_color
    module function particle_set_get_momenta_all (particle_set) result (p)
      class(particle_set_t), intent(in) :: particle_set
      type(vector4_t), dimension(:), allocatable :: p
    end function particle_set_get_momenta_all
    module function particle_set_get_momenta_indices (particle_set, indices) result (p)
       type(vector4_t), dimension(:), allocatable :: p
       class(particle_set_t), intent(in) :: particle_set
       integer, intent(in), dimension(:), allocatable :: indices
    end function particle_set_get_momenta_indices
    module subroutine particle_set_set_momentum_single &
         (particle_set, i, p, p2, on_shell)
      class(particle_set_t), intent(inout) :: particle_set
      integer, intent(in) :: i
      type(vector4_t), intent(in) :: p
      real(default), intent(in), optional :: p2
      logical, intent(in), optional :: on_shell
    end subroutine particle_set_set_momentum_single
    module subroutine particle_set_set_momentum_indices &
         (particle_set, indices, p, p2, on_shell)
      class(particle_set_t), intent(inout) :: particle_set
      integer, dimension(:), intent(in) :: indices
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), dimension(:), intent(in), optional :: p2
      logical, intent(in), optional :: on_shell
    end subroutine particle_set_set_momentum_indices
    module subroutine particle_set_set_momentum_all (particle_set, p, p2, on_shell)
      class(particle_set_t), intent(inout) :: particle_set
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), dimension(:), intent(in), optional :: p2
      logical, intent(in), optional :: on_shell
    end subroutine particle_set_set_momentum_all
    module subroutine particle_set_recover_momentum (particle_set, i)
      class(particle_set_t), intent(inout) :: particle_set
      integer, intent(in) :: i
    end subroutine particle_set_recover_momentum
    module subroutine particle_set_replace_incoming_momenta (particle_set, p)
      class(particle_set_t), intent(inout) :: particle_set
      type(vector4_t), intent(in), dimension(:) :: p
    end subroutine particle_set_replace_incoming_momenta
    module subroutine particle_set_replace_outgoing_momenta (particle_set, p)
      class(particle_set_t), intent(inout) :: particle_set
      type(vector4_t), intent(in), dimension(:) :: p
    end subroutine particle_set_replace_outgoing_momenta
    module function particle_set_get_outgoing_momenta (particle_set) result (p)
      class(particle_set_t), intent(in) :: particle_set
      type(vector4_t), dimension(:), allocatable :: p
    end function particle_set_get_outgoing_momenta
    module subroutine particle_set_parent_add_child (particle_set, parent, child)
      class(particle_set_t), intent(inout) :: particle_set
      integer, intent(in) :: parent, child
    end subroutine particle_set_parent_add_child
    module subroutine particle_set_build_radiation (particle_set, p_radiated, &
         emitter, flv_radiated, model, r_color)
      class(particle_set_t), intent(inout) :: particle_set
      type(vector4_t), intent(in), dimension(:) :: p_radiated
      integer, intent(in) :: emitter
      integer, intent(in), dimension(:) :: flv_radiated
      class(model_data_t), intent(in), target :: model
      real(default), intent(in) :: r_color
    end subroutine particle_set_build_radiation
    module subroutine particle_set_write &
      (particle_set, unit, testflag, summary, compressed)
      class(particle_set_t), intent(in) :: particle_set
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag, summary, compressed
    end subroutine particle_set_write
    module subroutine particle_set_write_raw (particle_set, u)
      class(particle_set_t), intent(in) :: particle_set
      integer, intent(in) :: u
    end subroutine particle_set_write_raw
    module subroutine particle_set_read_raw (particle_set, u, iostat)
      class(particle_set_t), intent(out) :: particle_set
      integer, intent(in) :: u
      integer, intent(out) :: iostat
    end subroutine particle_set_read_raw
    module function particle_set_get_real_parents (pset, i, keep_beams) result (parent)
      integer, dimension(:), allocatable :: parent
      class(particle_set_t), intent(in) :: pset
      integer, intent(in) :: i
      logical, intent(in), optional :: keep_beams
    end function particle_set_get_real_parents
    module function particle_set_get_real_children (pset, i, keep_beams) result (child)
      integer, dimension(:), allocatable :: child
      class(particle_set_t), intent(in) :: pset
      integer, intent(in) :: i
      logical, intent(in), optional :: keep_beams
    end function particle_set_get_real_children
    module function particle_set_get_n_beam (pset) result (n_beam)
       class(particle_set_t), intent(in) :: pset
       integer :: n_beam
    end function particle_set_get_n_beam
    module function particle_set_get_n_in (pset) result (n_in)
       class(particle_set_t), intent(in) :: pset
       integer :: n_in
    end function particle_set_get_n_in
    module function particle_set_get_n_vir (pset) result (n_vir)
       class(particle_set_t), intent(in) :: pset
       integer :: n_vir
     end function particle_set_get_n_vir
    module function particle_set_get_n_out (pset) result (n_out)
       class(particle_set_t), intent(in) :: pset
       integer :: n_out
    end function particle_set_get_n_out
    module function particle_set_get_n_tot (pset) result (n_tot)
       class(particle_set_t), intent(in) :: pset
       integer :: n_tot
    end function particle_set_get_n_tot
    module function particle_set_get_n_remnants (pset) result (n_remn)
      class(particle_set_t), intent(in) :: pset
      integer :: n_remn
    end function particle_set_get_n_remnants
    module function particle_set_get_particle (pset, index) result (particle)
      class(particle_set_t), intent(in) :: pset
      integer, intent(in) :: index
      type(particle_t) :: particle
    end function particle_set_get_particle
    pure module function particle_set_get_indices (pset, mask) result (finals)
      integer, dimension(:), allocatable :: finals
      class(particle_set_t), intent(in) :: pset
      logical, dimension(:), intent(in) :: mask
    end function particle_set_get_indices
    module function particle_set_get_in_and_out_momenta (pset) result (phs_point)
      type(phs_point_t) :: phs_point
      class(particle_set_t), intent(in) :: pset
    end function particle_set_get_in_and_out_momenta
    module subroutine particle_set_without_hadronic_remnants &
           (particle_set, particles, n_particles, n_extra)
      class(particle_set_t), intent(inout) :: particle_set
      type(particle_t), dimension(:), allocatable, intent(out) :: particles
      integer, intent(out) :: n_particles
      integer, intent(in) :: n_extra
    end subroutine particle_set_without_hadronic_remnants
    module subroutine particle_set_without_remnants &
           (particle_set, particles, n_particles, n_extra)
      class(particle_set_t), intent(inout) :: particle_set
      type(particle_t), dimension(:), allocatable, intent(out) :: particles
      integer, intent(in) :: n_extra
      integer, intent(out) :: n_particles
    end subroutine particle_set_without_remnants
    pure module function particle_set_find_particle (particle_set, pdg, &
         momentum, abs_smallness, rel_smallness) result (idx)
      integer :: idx
      class(particle_set_t), intent(in) :: particle_set
      integer, intent(in) :: pdg
      type(vector4_t), intent(in) :: momentum
      real(default), intent(in), optional :: abs_smallness, rel_smallness
    end function particle_set_find_particle
    pure module function particle_set_reverse_find_particle &
         (particle_set, pdg, momentum, abs_smallness, rel_smallness) result (idx)
      integer :: idx
      class(particle_set_t), intent(in) :: particle_set
      integer, intent(in) :: pdg
      type(vector4_t), intent(in) :: momentum
      real(default), intent(in), optional :: abs_smallness, rel_smallness
    end function particle_set_reverse_find_particle
    module subroutine particle_set_remove_duplicates (particle_set, smallness)
      class(particle_set_t), intent(inout) :: particle_set
      real(default), intent(in) :: smallness
    end subroutine particle_set_remove_duplicates
    module subroutine particle_set_reset_status (particle_set, index, status)
      class(particle_set_t), intent(inout) :: particle_set
      integer, dimension(:), intent(in) :: index
      integer, intent(in) :: status
    end subroutine particle_set_reset_status
    module subroutine particle_set_reduce (pset_in, pset_out, keep_beams)
      class(particle_set_t), intent(in) :: pset_in
      type(particle_set_t), intent(out) :: pset_out
      logical, intent(in), optional :: keep_beams
    end subroutine particle_set_reduce
    module subroutine particle_set_filter_particles &
         (pset_in, pset_out, keep_beams, real_parents, keep_virtuals)
      class(particle_set_t), intent(in) :: pset_in
      type(particle_set_t), intent(out) :: pset_out
      logical, intent(in), optional :: keep_beams, real_parents, keep_virtuals
    end subroutine particle_set_filter_particles
    module subroutine particle_set_to_hepevt_form (pset_in, pset_out)
      class(particle_set_t), intent(in) :: pset_in
      type(particle_set_t), intent(out) :: pset_out
    end subroutine particle_set_to_hepevt_form
    module subroutine particle_set_fill_interaction &
         (pset, int, n_in, recover_beams, check_match, state_flv, success)
      class(particle_set_t), intent(in) :: pset
      type(interaction_t), intent(inout) :: int
      integer, intent(in) :: n_in
      logical, intent(in), optional :: recover_beams, check_match
      type(state_flv_content_t), intent(in), optional :: state_flv
      logical, intent(out), optional :: success
    end subroutine particle_set_fill_interaction
    module subroutine particle_set_assign_vertices &
         (particle_set, v_from, v_to, n_vertices)
      class(particle_set_t), intent(in) :: particle_set
      integer, dimension(:), intent(out) :: v_from, v_to
      integer, intent(out) :: n_vertices
    end subroutine particle_set_assign_vertices
    module subroutine particle_set_to_subevt (particle_set, subevt, colorize)
      class(particle_set_t), intent(in) :: particle_set
      type(subevt_t), intent(out) :: subevt
      logical, intent(in), optional :: colorize
    end subroutine particle_set_to_subevt
    module subroutine particle_set_replace (particle_set, newprt)
      class(particle_set_t), intent(inout) :: particle_set
      type(particle_t), intent(in), dimension(:), allocatable :: newprt
    end subroutine particle_set_replace
    module subroutine particle_set_order_color_lines (pset_out, pset_in)
      class(particle_set_t), intent(inout) :: pset_out
      type(particle_set_t), intent(in) :: pset_in
    end subroutine particle_set_order_color_lines
    module subroutine pacify_particle (prt)
      class(particle_t), intent(inout) :: prt
    end subroutine pacify_particle
    module subroutine pacify_particle_set (pset)
      class(particle_set_t), intent(inout) :: pset
    end subroutine pacify_particle_set
  end interface

end module particles
