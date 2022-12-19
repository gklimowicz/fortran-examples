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

module flavors

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use physics_defs, only: UNDEFINED
  use model_data
  use colors, only: color_t

  implicit none
  private

  public :: flavor_t
  public :: flavor_write_array
  public :: operator(.merge.)
  public :: color_from_flavor

  type :: flavor_t
     private
     integer :: f = UNDEFINED
     logical :: hard_process = .false.
     logical :: radiated = .false.
     type(field_data_t), pointer :: field_data => null ()
   contains
     generic :: init => &
          flavor_init_empty, &
          flavor_init, &
          flavor_init_field_data, &
          flavor_init_model, &
          flavor_init_model_alt, &
          flavor_init_name_model
     procedure, private :: flavor_init_empty
     procedure, private :: flavor_init
     procedure, private :: flavor_init_field_data
     procedure, private :: flavor_init_model
     procedure, private :: flavor_init_model_alt
     procedure, private :: flavor_init_name_model
     procedure :: tag_radiated => flavor_tag_radiated
     procedure :: tag_hard_process => flavor_tag_hard_process
     procedure :: undefine => flavor_undefine
     procedure :: write => flavor_write
     procedure :: write_raw => flavor_write_raw
     procedure :: read_raw => flavor_read_raw
     procedure :: set_model => flavor_set_model_single
     procedure :: is_defined => flavor_is_defined
     procedure :: is_valid => flavor_is_valid
     procedure :: is_associated => flavor_is_associated
     procedure :: is_radiated => flavor_is_radiated
     procedure :: is_hard_process => flavor_is_hard_process
     procedure :: get_pdg => flavor_get_pdg
     procedure :: get_pdg_anti => flavor_get_pdg_anti
     procedure :: get_pdg_abs => flavor_get_pdg_abs
     procedure :: is_visible => flavor_is_visible
     procedure :: is_parton => flavor_is_parton
     procedure :: is_beam_remnant => flavor_is_beam_remnant
     procedure :: is_gauge => flavor_is_gauge
     procedure :: is_left_handed => flavor_is_left_handed
     procedure :: is_right_handed => flavor_is_right_handed
     procedure :: is_antiparticle => flavor_is_antiparticle
     procedure :: has_antiparticle => flavor_has_antiparticle
     procedure :: is_stable => flavor_is_stable
     procedure :: get_decays => flavor_get_decays
     procedure :: decays_isotropically => flavor_decays_isotropically
     procedure :: decays_diagonal => flavor_decays_diagonal
     procedure :: has_decay_helicity => flavor_has_decay_helicity
     procedure :: get_decay_helicity => flavor_get_decay_helicity
     procedure :: is_polarized => flavor_is_polarized
     procedure :: get_name => flavor_get_name
     procedure :: get_tex_name => flavor_get_tex_name
     procedure :: get_spin_type => flavor_get_spin_type
     procedure :: get_multiplicity => flavor_get_multiplicity
     procedure :: get_isospin_type => flavor_get_isospin_type
     procedure :: get_charge_type => flavor_get_charge_type
     procedure :: get_color_type => flavor_get_color_type
     procedure :: get_charge => flavor_get_charge
     procedure :: get_mass => flavor_get_mass
     procedure :: get_width => flavor_get_width
     procedure :: get_isospin => flavor_get_isospin
     generic :: operator(.match.) => flavor_match
     generic :: operator(==) => flavor_eq
     generic :: operator(/=) => flavor_neq
     procedure, private :: flavor_match
     procedure, private :: flavor_eq
     procedure, private :: flavor_neq
     procedure :: anti => flavor_anti
  end type flavor_t


  interface operator(.merge.)
     module procedure merge_flavors0
     module procedure merge_flavors1
  end interface

  interface color_from_flavor
     module procedure color_from_flavor0
     module procedure color_from_flavor1
  end interface

  interface
    elemental module subroutine flavor_init_empty (flv)
      class(flavor_t), intent(inout) :: flv
    end subroutine flavor_init_empty
    elemental module subroutine flavor_init (flv, f)
      class(flavor_t), intent(inout) :: flv
      integer, intent(in) :: f
    end subroutine flavor_init
    impure elemental module subroutine flavor_init_field_data (flv, field_data)
      class(flavor_t), intent(inout) :: flv
      type(field_data_t), intent(in), target :: field_data
    end subroutine flavor_init_field_data
    impure elemental module subroutine flavor_init_model (flv, f, model)
      class(flavor_t), intent(inout) :: flv
      integer, intent(in) :: f
      class(model_data_t), intent(in), target :: model
    end subroutine flavor_init_model
    impure elemental module subroutine flavor_init_model_alt (flv, f, model, alt_model)
      class(flavor_t), intent(inout) :: flv
      integer, intent(in) :: f
      class(model_data_t), intent(in), target :: model, alt_model
    end subroutine flavor_init_model_alt
    impure elemental module subroutine flavor_init_name_model (flv, name, model)
      class(flavor_t), intent(inout) :: flv
      type(string_t), intent(in) :: name
      class(model_data_t), intent(in), target :: model
    end subroutine flavor_init_name_model
    elemental module subroutine flavor_tag_radiated (flv)
      class(flavor_t), intent(inout) :: flv
    end subroutine flavor_tag_radiated
    elemental module subroutine flavor_tag_hard_process (flv, hard)
      class(flavor_t), intent(inout) :: flv
      logical, intent(in), optional :: hard
    end subroutine flavor_tag_hard_process
    elemental module subroutine flavor_undefine (flv)
      class(flavor_t), intent(inout) :: flv
    end subroutine flavor_undefine
    module subroutine flavor_write (flv, unit)
      class(flavor_t), intent(in) :: flv
      integer, intent(in), optional :: unit
    end subroutine flavor_write
    module subroutine flavor_write_array (flv, unit)
      type(flavor_t), intent(in), dimension(:) :: flv
      integer, intent(in), optional :: unit
    end subroutine flavor_write_array
    module subroutine flavor_write_raw (flv, u)
      class(flavor_t), intent(in) :: flv
      integer, intent(in) :: u
    end subroutine flavor_write_raw
    module subroutine flavor_read_raw (flv, u, iostat)
      class(flavor_t), intent(out) :: flv
      integer, intent(in) :: u
      integer, intent(out), optional :: iostat
    end subroutine flavor_read_raw
    impure elemental module subroutine flavor_set_model_single (flv, model)
      class(flavor_t), intent(inout) :: flv
      class(model_data_t), intent(in), target :: model
    end subroutine flavor_set_model_single
    elemental module function flavor_is_defined (flv) result (defined)
      class(flavor_t), intent(in) :: flv
      logical :: defined
    end function flavor_is_defined
    elemental module function flavor_is_valid (flv) result (valid)
      class(flavor_t), intent(in) :: flv
      logical :: valid
    end function flavor_is_valid
    elemental module function flavor_is_associated (flv) result (flag)
      class(flavor_t), intent(in) :: flv
      logical :: flag
    end function flavor_is_associated
    elemental module function flavor_is_radiated (flv) result (flag)
      class(flavor_t), intent(in) :: flv
      logical :: flag
    end function flavor_is_radiated
    elemental module function flavor_is_hard_process (flv) result (flag)
      class(flavor_t), intent(in) :: flv
      logical :: flag
    end function flavor_is_hard_process
    elemental module function flavor_get_pdg (flv) result (f)
      integer :: f
      class(flavor_t), intent(in) :: flv
    end function flavor_get_pdg
    elemental module function flavor_get_pdg_anti (flv) result (f)
      integer :: f
      class(flavor_t), intent(in) :: flv
    end function flavor_get_pdg_anti
    elemental module function flavor_get_pdg_abs (flv) result (f)
      integer :: f
      class(flavor_t), intent(in) :: flv
    end function flavor_get_pdg_abs
    elemental module function flavor_is_visible (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_is_visible
    elemental module function flavor_is_parton (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_is_parton
    elemental module function flavor_is_beam_remnant (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_is_beam_remnant
    elemental module function flavor_is_gauge (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_is_gauge
    elemental module function flavor_is_left_handed (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_is_left_handed
    elemental module function flavor_is_right_handed (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_is_right_handed
    elemental module function flavor_is_antiparticle (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_is_antiparticle
    elemental module function flavor_has_antiparticle (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_has_antiparticle
    elemental module function flavor_is_stable (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_is_stable
    module subroutine flavor_get_decays (flv, decay)
      class(flavor_t), intent(in) :: flv
      type(string_t), dimension(:), intent(out), allocatable :: decay
      logical :: anti
    end subroutine flavor_get_decays
    elemental module function flavor_decays_isotropically (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_decays_isotropically
    elemental module function flavor_decays_diagonal (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_decays_diagonal
    elemental module function flavor_has_decay_helicity (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_has_decay_helicity
    elemental module function flavor_get_decay_helicity (flv) result (hel)
      integer :: hel
      class(flavor_t), intent(in) :: flv
    end function flavor_get_decay_helicity
    elemental module function flavor_is_polarized (flv) result (flag)
      logical :: flag
      class(flavor_t), intent(in) :: flv
    end function flavor_is_polarized
    elemental module function flavor_get_name (flv) result (name)
      type(string_t) :: name
      class(flavor_t), intent(in) :: flv
    end function flavor_get_name
    elemental module function flavor_get_tex_name (flv) result (name)
      type(string_t) :: name
      class(flavor_t), intent(in) :: flv
    end function flavor_get_tex_name
    elemental module function flavor_get_spin_type (flv) result (type)
      integer :: type
      class(flavor_t), intent(in) :: flv
    end function flavor_get_spin_type
    elemental module function flavor_get_multiplicity (flv) result (type)
      integer :: type
      class(flavor_t), intent(in) :: flv
    end function flavor_get_multiplicity
    elemental module function flavor_get_isospin_type (flv) result (type)
      integer :: type
      class(flavor_t), intent(in) :: flv
    end function flavor_get_isospin_type
    elemental module function flavor_get_charge_type (flv) result (type)
      integer :: type
      class(flavor_t), intent(in) :: flv
    end function flavor_get_charge_type
    elemental module function flavor_get_color_type (flv) result (type)
      integer :: type
      class(flavor_t), intent(in) :: flv
    end function flavor_get_color_type
    elemental module function flavor_get_charge (flv) result (charge)
      real(default) :: charge
      class(flavor_t), intent(in) :: flv
    end function flavor_get_charge
    elemental module function flavor_get_mass (flv) result (mass)
      real(default) :: mass
      class(flavor_t), intent(in) :: flv
    end function flavor_get_mass
    elemental module function flavor_get_width (flv) result (width)
      real(default) :: width
      class(flavor_t), intent(in) :: flv
    end function flavor_get_width
    elemental module function flavor_get_isospin (flv) result (isospin)
      real(default) :: isospin
      class(flavor_t), intent(in) :: flv
    end function flavor_get_isospin
    elemental module function flavor_match (flv1, flv2) result (eq)
      logical :: eq
      class(flavor_t), intent(in) :: flv1, flv2
    end function flavor_match
    elemental module function flavor_eq (flv1, flv2) result (eq)
      logical :: eq
      class(flavor_t), intent(in) :: flv1, flv2
    end function flavor_eq
    elemental module function flavor_neq (flv1, flv2) result (neq)
      logical :: neq
      class(flavor_t), intent(in) :: flv1, flv2
    end function flavor_neq
    module function merge_flavors0 (flv1, flv2) result (flv)
      type(flavor_t) :: flv
      type(flavor_t), intent(in) :: flv1, flv2
    end function merge_flavors0
    module function merge_flavors1 (flv1, flv2) result (flv)
      type(flavor_t), dimension(:), intent(in) :: flv1, flv2
      type(flavor_t), dimension(size(flv1)) :: flv
    end function merge_flavors1
    module function color_from_flavor0 (flv, c_seed, reverse) result (col)
      type(color_t) :: col
      type(flavor_t), intent(in) :: flv
      integer, intent(in), optional :: c_seed
      logical, intent(in), optional :: reverse
    end function color_from_flavor0
    module function color_from_flavor1 (flv, c_seed, reverse) result (col)
      type(flavor_t), dimension(:), intent(in) :: flv
      integer, intent(in), optional :: c_seed
      logical, intent(in), optional :: reverse
      type(color_t), dimension(size(flv)) :: col
    end function color_from_flavor1
    module function flavor_anti (flv) result (aflv)
      type(flavor_t) :: aflv
      class(flavor_t), intent(in) :: flv
    end function flavor_anti
  end interface

end module flavors
