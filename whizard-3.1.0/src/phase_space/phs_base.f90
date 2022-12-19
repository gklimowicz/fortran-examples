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

module phs_base

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use model_data
  use flavors
  use process_constants

  implicit none
  private

  public :: channel_prop_t
  public :: resonance_t
  public :: on_shell_t
  public :: phs_channel_t
  public :: phs_channel_collection_t
  public :: phs_config_t
  public :: phs_t
  public :: compute_kinematics_solid_angle
  public :: inverse_kinematics_solid_angle
  public :: pacify

  type, abstract :: channel_prop_t
   contains
     procedure (channel_prop_to_string), deferred :: to_string
     generic :: operator (==) => is_equal
     procedure (channel_eq), deferred :: is_equal
  end type channel_prop_t

  type, extends (channel_prop_t) :: resonance_t
     real(default) :: mass = 0
     real(default) :: width = 0
   contains
     procedure :: to_string => resonance_to_string
     procedure :: is_equal => resonance_is_equal
  end type resonance_t

  type, extends (channel_prop_t) :: on_shell_t
     real(default) :: mass = 0
   contains
     procedure :: to_string => on_shell_to_string
     procedure :: is_equal => on_shell_is_equal
  end type on_shell_t

  type :: phs_equivalence_t
     integer :: c = 0
     integer, dimension(:), allocatable :: perm
     integer, dimension(:), allocatable :: mode
   contains
     procedure :: write => phs_equivalence_write
     procedure :: init => phs_equivalence_init
  end type phs_equivalence_t

  integer, parameter, public :: &
       EQ_IDENTITY = 0, EQ_INVERT = 1, EQ_SYMMETRIC = 2, EQ_INVARIANT = 3

  character, dimension(0:3), parameter :: TAG = ["+", "-", ":", "x"]

  type :: phs_channel_t
     class(channel_prop_t), allocatable :: prop
     integer :: sf_channel = 1
     type(phs_equivalence_t), dimension(:), allocatable :: eq
   contains
     procedure :: write => phs_channel_write
     procedure :: set_resonant => channel_set_resonant
     procedure :: set_on_shell => channel_set_on_shell
  end type phs_channel_t

  type :: prop_entry_t
     integer :: i = 0
     class(channel_prop_t), allocatable :: prop
     type(prop_entry_t), pointer :: next => null ()
  end type prop_entry_t

  type :: phs_channel_collection_t
     integer :: n = 0
     type(prop_entry_t), pointer :: first => null ()
   contains
     procedure :: final => phs_channel_collection_final
     procedure :: write => phs_channel_collection_write
     procedure :: push => phs_channel_collection_push
     procedure :: get_n => phs_channel_collection_get_n
     procedure :: get_entry => phs_channel_collection_get_entry
  end type phs_channel_collection_t

  type, abstract :: phs_config_t
     ! private
     type(string_t) :: id
     integer :: n_in = 0
     integer :: n_out = 0
     integer :: n_tot = 0
     integer :: n_state = 0
     integer :: n_par = 0
     integer :: n_channel = 0
     real(default) :: sqrts = 0
     logical :: sqrts_fixed = .true.
     logical :: lab_is_cm = .true.
     logical :: azimuthal_dependence = .false.
     integer, dimension(:), allocatable :: dim_flat
     logical :: provides_equivalences = .false.
     logical :: provides_chains = .false.
     logical :: vis_channels = .false.
     integer, dimension(:), allocatable :: chain
     class(model_data_t), pointer :: model => null ()
     type(flavor_t), dimension(:,:), allocatable :: flv
     type(phs_channel_t), dimension(:), allocatable :: channel
     character(32) :: md5sum_process = ""
     character(32) :: md5sum_model_par = ""
     character(32) :: md5sum_phs_config = ""
     integer :: nlo_type
   contains
     procedure (phs_config_final), deferred :: final
     procedure (phs_config_write), deferred :: write
     procedure :: base_write => phs_config_write
     procedure :: init => phs_config_init
     procedure (phs_config_configure), deferred :: configure
     procedure :: set_sf_channel => phs_config_set_sf_channel
     procedure :: collect_channels => phs_config_collect_channels
     procedure :: compute_md5sum => phs_config_compute_md5sum
     procedure (phs_startup_message), deferred :: startup_message
     procedure :: base_startup_message => phs_startup_message
     procedure (phs_config_allocate_instance), nopass, deferred :: &
          allocate_instance
     procedure :: get_n_par => phs_config_get_n_par
     procedure :: get_flat_dimensions => phs_config_get_flat_dimensions
     procedure :: get_n_channel => phs_config_get_n_channel
     procedure :: get_sf_channel => phs_config_get_sf_channel
     procedure :: get_masses_in => phs_config_get_masses_in
     procedure :: get_md5sum => phs_config_get_md5sum
  end type phs_config_t

  type, abstract :: phs_t
     class(phs_config_t), pointer :: config => null ()
     logical :: r_defined = .false.
     integer :: selected_channel = 0
     logical, dimension(:), allocatable :: active_channel
     real(default), dimension(:,:), allocatable :: r
     real(default), dimension(:), allocatable :: f
     real(default), dimension(:), allocatable :: m_in
     real(default), dimension(:), allocatable :: m_out
     real(default) :: flux = 0
     real(default) :: volume = 0
     type(lorentz_transformation_t) :: lt_cm_to_lab
     logical :: p_defined = .false.
     real(default) :: sqrts_hat = 0
     type(vector4_t), dimension(:), allocatable :: p
     logical :: q_defined = .false.
     type(vector4_t), dimension(:), allocatable :: q
   contains
     procedure (phs_write), deferred :: write
     procedure :: base_write => phs_base_write
     procedure (phs_final), deferred :: final
     procedure (phs_init), deferred :: init
     procedure :: base_init => phs_base_init
     procedure :: select_channel => phs_base_select_channel
     procedure :: set_incoming_momenta => phs_set_incoming_momenta
     procedure :: set_outgoing_momenta => phs_set_outgoing_momenta
     procedure :: get_outgoing_momenta => phs_get_outgoing_momenta
     procedure :: lab_is_cm => phs_lab_is_cm
     procedure :: get_n_tot => phs_get_n_tot
     procedure :: set_lorentz_transformation => phs_set_lorentz_transformation
     procedure :: get_lorentz_transformation => phs_get_lorentz_transformation
     procedure :: get_mcpar => phs_get_mcpar
     procedure :: get_f => phs_get_f
     procedure :: get_overall_factor => phs_get_overall_factor
     procedure :: compute_flux => phs_compute_flux
     procedure :: compute_base_flux => phs_compute_flux
     procedure (phs_evaluate_selected_channel), deferred :: &
          evaluate_selected_channel
     procedure (phs_evaluate_other_channels), deferred :: &
          evaluate_other_channels
     procedure (phs_inverse), deferred :: inverse
     procedure :: get_sqrts => phs_get_sqrts
  end type phs_t


  abstract interface
     function channel_prop_to_string (object) result (string)
       import
       class(channel_prop_t), intent(in) :: object
       type(string_t) :: string
     end function channel_prop_to_string
  end interface

  abstract interface
     function channel_eq (prop1, prop2) result (flag)
       import
       class(channel_prop_t), intent(in) :: prop1, prop2
       logical :: flag
     end function channel_eq
  end interface

  abstract interface
     subroutine phs_config_final (object)
       import
       class(phs_config_t), intent(inout) :: object
     end subroutine phs_config_final
  end interface

  abstract interface
     subroutine phs_config_configure (phs_config, sqrts, &
          sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, ignore_mismatch, &
          nlo_type, subdir)
       import
       class(phs_config_t), intent(inout) :: phs_config
       real(default), intent(in) :: sqrts
       logical, intent(in), optional :: sqrts_fixed
       logical, intent(in), optional :: lab_is_cm
       logical, intent(in), optional :: azimuthal_dependence
       logical, intent(in), optional :: rebuild
       logical, intent(in), optional :: ignore_mismatch
       integer, intent(in), optional :: nlo_type
       type(string_t), intent(in), optional :: subdir
     end subroutine phs_config_configure
  end interface

  abstract interface
     subroutine phs_config_allocate_instance (phs)
       import
       class(phs_t), intent(inout), pointer :: phs
     end subroutine phs_config_allocate_instance
  end interface

  abstract interface
     subroutine phs_write (object, unit, verbose)
       import
       class(phs_t), intent(in) :: object
       integer, intent(in), optional :: unit
       logical, intent(in), optional :: verbose
     end subroutine phs_write
  end interface

  abstract interface
     subroutine phs_final (object)
       import
       class(phs_t), intent(inout) :: object
     end subroutine phs_final
  end interface

  abstract interface
     subroutine phs_init (phs, phs_config)
       import
       class(phs_t), intent(out) :: phs
       class(phs_config_t), intent(in), target :: phs_config
     end subroutine phs_init
  end interface

  abstract interface
     subroutine phs_evaluate_selected_channel (phs, c_in, r_in)
       import
       class(phs_t), intent(inout) :: phs
       integer, intent(in) :: c_in
       real(default), dimension(:), intent(in) :: r_in
     end subroutine phs_evaluate_selected_channel
  end interface

  abstract interface
     subroutine phs_evaluate_other_channels (phs, c_in)
       import
       class(phs_t), intent(inout) :: phs
       integer, intent(in) :: c_in
     end subroutine phs_evaluate_other_channels
  end interface

  abstract interface
     subroutine phs_inverse (phs)
       import
       class(phs_t), intent(inout) :: phs
     end subroutine phs_inverse
  end interface

  interface pacify
     module procedure pacify_phs
  end interface pacify


  interface
    module function resonance_to_string (object) result (string)
      class(resonance_t), intent(in) :: object
      type(string_t) :: string
    end function resonance_to_string
    module function resonance_is_equal (prop1, prop2) result (flag)
      class(resonance_t), intent(in) :: prop1
      class(channel_prop_t), intent(in) :: prop2
      logical :: flag
    end function resonance_is_equal
    module function on_shell_to_string (object) result (string)
      class(on_shell_t), intent(in) :: object
      type(string_t) :: string
    end function on_shell_to_string
    module function on_shell_is_equal (prop1, prop2) result (flag)
      class(on_shell_t), intent(in) :: prop1
      class(channel_prop_t), intent(in) :: prop2
      logical :: flag
    end function on_shell_is_equal
    module subroutine phs_equivalence_write (object, unit)
      class(phs_equivalence_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine phs_equivalence_write
    module subroutine phs_equivalence_init (eq, n_dim)
      class(phs_equivalence_t), intent(out) :: eq
      integer, intent(in) :: n_dim
    end subroutine phs_equivalence_init
    module subroutine phs_channel_write (object, unit)
      class(phs_channel_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine phs_channel_write
    module subroutine phs_channel_collection_final (object)
      class(phs_channel_collection_t), intent(inout) :: object
    end subroutine phs_channel_collection_final
    module subroutine phs_channel_collection_write (object, unit)
      class(phs_channel_collection_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine phs_channel_collection_write
    module subroutine phs_channel_collection_push (coll, channel)
      class(phs_channel_collection_t), intent(inout) :: coll
      type(phs_channel_t), intent(inout) :: channel
    end subroutine phs_channel_collection_push
    module function phs_channel_collection_get_n (coll) result (n)
      class(phs_channel_collection_t), intent(in) :: coll
      integer :: n
    end function phs_channel_collection_get_n
    module subroutine phs_channel_collection_get_entry (coll, i, prop)
      class(phs_channel_collection_t), intent(in) :: coll
      integer, intent(in) :: i
      class(channel_prop_t), intent(out), allocatable :: prop
    end subroutine phs_channel_collection_get_entry
    module subroutine phs_config_write (object, unit, include_id)
      class(phs_config_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: include_id
    end subroutine phs_config_write
    module subroutine phs_config_init (phs_config, data, model)
      class(phs_config_t), intent(inout) :: phs_config
      type(process_constants_t), intent(in) :: data
      class(model_data_t), intent(in), target :: model
    end subroutine phs_config_init
    module subroutine phs_config_set_sf_channel (phs_config, sf_channel)
      class(phs_config_t), intent(inout) :: phs_config
      integer, dimension(:), intent(in) :: sf_channel
    end subroutine phs_config_set_sf_channel
    module subroutine phs_config_collect_channels (phs_config, coll)
      class(phs_config_t), intent(inout) :: phs_config
      type(phs_channel_collection_t), intent(inout) :: coll
    end subroutine phs_config_collect_channels
    module subroutine phs_config_compute_md5sum (phs_config, include_id)
      class(phs_config_t), intent(inout) :: phs_config
      logical, intent(in), optional :: include_id
    end subroutine phs_config_compute_md5sum
    module subroutine phs_startup_message (phs_config, unit)
      class(phs_config_t), intent(in) :: phs_config
      integer, intent(in), optional :: unit
    end subroutine phs_startup_message
    module function phs_config_get_n_par (phs_config) result (n)
      class(phs_config_t), intent(in) :: phs_config
      integer :: n
    end function phs_config_get_n_par
    module function phs_config_get_flat_dimensions &
         (phs_config) result (dim_flat)
      class(phs_config_t), intent(in) :: phs_config
      integer, dimension(:), allocatable :: dim_flat
    end function phs_config_get_flat_dimensions
    module function phs_config_get_n_channel (phs_config) result (n)
      class(phs_config_t), intent(in) :: phs_config
      integer :: n
    end function phs_config_get_n_channel
    module function phs_config_get_sf_channel (phs_config, c) result (c_sf)
      class(phs_config_t), intent(in) :: phs_config
      integer, intent(in) :: c
      integer :: c_sf
    end function phs_config_get_sf_channel
    module subroutine phs_config_get_masses_in (phs_config, m)
      class(phs_config_t), intent(in) :: phs_config
      real(default), dimension(:), intent(out) :: m
    end subroutine phs_config_get_masses_in
    module function phs_config_get_md5sum (phs_config) result (md5sum)
      class(phs_config_t), intent(in) :: phs_config
      character(32) :: md5sum
    end function phs_config_get_md5sum
    module subroutine phs_base_write (object, unit)
      class(phs_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine phs_base_write
    module subroutine phs_base_init (phs, phs_config)
      class(phs_t), intent(out) :: phs
      class(phs_config_t), intent(in), target :: phs_config
    end subroutine phs_base_init
    module subroutine phs_base_select_channel (phs, channel)
      class(phs_t), intent(inout) :: phs
      integer, intent(in), optional :: channel
    end subroutine phs_base_select_channel
    module subroutine phs_set_incoming_momenta (phs, p)
      class(phs_t), intent(inout) :: phs
      type(vector4_t), dimension(:), intent(in) :: p
    end subroutine phs_set_incoming_momenta
    module subroutine phs_set_outgoing_momenta (phs, q)
      class(phs_t), intent(inout) :: phs
      type(vector4_t), dimension(:), intent(in) :: q
    end subroutine phs_set_outgoing_momenta
    module subroutine phs_get_outgoing_momenta (phs, q)
      class(phs_t), intent(in) :: phs
      type(vector4_t), dimension(:), intent(out) :: q
    end subroutine phs_get_outgoing_momenta
    module function phs_lab_is_cm (phs) result (lab_is_cm)
      logical :: lab_is_cm
      class(phs_t), intent(in) :: phs
    end function phs_lab_is_cm
    elemental module function phs_get_n_tot (phs) result (n_tot)
      integer :: n_tot
      class(phs_t), intent(in) :: phs
    end function phs_get_n_tot
    module subroutine phs_set_lorentz_transformation (phs, lt)
      class(phs_t), intent(inout) :: phs
      type(lorentz_transformation_t), intent(in) :: lt
    end subroutine phs_set_lorentz_transformation
    module function phs_get_lorentz_transformation (phs) result (lt)
      type(lorentz_transformation_t) :: lt
      class(phs_t), intent(in) :: phs
    end function phs_get_lorentz_transformation
    module subroutine phs_get_mcpar (phs, c, r)
      class(phs_t), intent(in) :: phs
      integer, intent(in) :: c
      real(default), dimension(:), intent(out) :: r
    end subroutine phs_get_mcpar
    module function phs_get_f (phs, c) result (f)
      class(phs_t), intent(in) :: phs
      integer, intent(in) :: c
      real(default) :: f
    end function phs_get_f
    module function phs_get_overall_factor (phs) result (f)
      class(phs_t), intent(in) :: phs
      real(default) :: f
    end function phs_get_overall_factor
    module subroutine phs_compute_flux (phs)
      class(phs_t), intent(inout) :: phs
    end subroutine phs_compute_flux
    module function phs_get_sqrts (phs) result (sqrts)
      real(default) :: sqrts
      class(phs_t), intent(in) :: phs
    end function phs_get_sqrts
    module subroutine compute_kinematics_solid_angle (p, q, x)
      type(vector4_t), dimension(2), intent(in) :: p
      type(vector4_t), dimension(2), intent(out) :: q
      real(default), dimension(2), intent(in) :: x
    end subroutine compute_kinematics_solid_angle
    module subroutine inverse_kinematics_solid_angle (p, q, x)
      type(vector4_t), dimension(:), intent(in) :: p
      type(vector4_t), dimension(2), intent(in) :: q
      real(default), dimension(2), intent(out) :: x
    end subroutine inverse_kinematics_solid_angle
    module subroutine pacify_phs (phs)
      class(phs_t), intent(inout) :: phs
    end subroutine pacify_phs
  end interface

contains

  subroutine channel_set_resonant (channel, mass, width)
    class(phs_channel_t), intent(inout) :: channel
    real(default), intent(in) :: mass, width
    allocate (resonance_t :: channel%prop)
    select type (prop => channel%prop)
    type is (resonance_t)
       prop%mass = mass
       prop%width = width
    end select
  end subroutine channel_set_resonant

  subroutine channel_set_on_shell (channel, mass)
    class(phs_channel_t), intent(inout) :: channel
    real(default), intent(in) :: mass
    allocate (on_shell_t :: channel%prop)
    select type (prop => channel%prop)
    type is (on_shell_t)
       prop%mass = mass
    end select
  end subroutine channel_set_on_shell


end module phs_base
