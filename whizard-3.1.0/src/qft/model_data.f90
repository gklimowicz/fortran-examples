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

module model_data

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds, only: default
  use kinds, only: i8, i32
  use kinds, only: c_default_float
  use iso_varying_string, string_t => varying_string
  use physics_defs, only: UNDEFINED, SCALAR

  implicit none
  private

  public :: modelpar_data_t
  public :: field_data_t
  public :: find_model
  public :: vertex_iterator_t
  public :: model_data_t

  type, abstract :: modelpar_data_t
     private
     type(string_t) :: name
   contains
     procedure :: write => par_write
     procedure :: show => par_show
     generic :: init => modelpar_data_init_real, modelpar_data_init_complex
     procedure, private :: modelpar_data_init_real
     procedure, private :: modelpar_data_init_complex
     generic :: assignment(=) => modelpar_data_set_real, modelpar_data_set_complex
     procedure, private :: modelpar_data_set_real
     procedure, private :: modelpar_data_set_complex
     procedure :: get_name => modelpar_data_get_name
     procedure, pass :: get_real => modelpar_data_get_real
     procedure, pass :: get_complex => modelpar_data_get_complex
     procedure :: get_real_ptr => modelpar_data_get_real_ptr
     procedure :: get_complex_ptr => modelpar_data_get_complex_ptr
  end type modelpar_data_t

  type, extends (modelpar_data_t) :: modelpar_real_t
     private
     real(default) :: value
  end type modelpar_real_t

  type, extends (modelpar_data_t) :: modelpar_complex_t
     private
     complex(default) :: value
  end type modelpar_complex_t

  type :: field_data_t
     private
     type(string_t) :: longname
     integer :: pdg = UNDEFINED
     logical :: visible = .true.
     logical :: parton = .false.
     logical :: gauge = .false.
     logical :: left_handed = .false.
     logical :: right_handed = .false.
     logical :: has_anti = .false.
     logical :: p_is_stable = .true.
     logical :: p_decays_isotropically = .false.
     logical :: p_decays_diagonal = .false.
     logical :: p_has_decay_helicity = .false.
     integer :: p_decay_helicity = 0
     logical :: a_is_stable = .true.
     logical :: a_decays_isotropically = .false.
     logical :: a_decays_diagonal = .false.
     logical :: a_has_decay_helicity = .false.
     integer :: a_decay_helicity = 0
     logical :: p_polarized = .false.
     logical :: a_polarized = .false.
     type(string_t), dimension(:), allocatable :: name, anti
     type(string_t) :: tex_name, tex_anti
     integer :: spin_type = UNDEFINED
     integer :: isospin_type = 1
     integer :: charge_type = 1
     integer :: color_type = 1
     real(default), pointer :: mass_val => null ()
     class(modelpar_data_t), pointer :: mass_data => null ()
     real(default), pointer :: width_val => null ()
     class(modelpar_data_t), pointer :: width_data => null ()
     integer :: multiplicity = 1
     type(string_t), dimension(:), allocatable :: p_decay
     type(string_t), dimension(:), allocatable :: a_decay
   contains
     procedure :: init => field_data_init
     procedure :: copy_from => field_data_copy_from
     procedure :: set => field_data_set
     procedure, private :: &
          set_multiplicity => field_data_set_multiplicity
     procedure, private :: set_mass => field_data_set_mass
     procedure, private :: set_width => field_data_set_width
     procedure :: freeze => field_data_freeze
     procedure :: write => field_data_write
     procedure :: write_decays => field_data_write_decays
     procedure :: show => field_data_show
     procedure :: get_pdg => field_data_get_pdg
     procedure :: get_pdg_anti => field_data_get_pdg_anti
     procedure :: is_visible => field_data_is_visible
     procedure :: is_parton => field_data_is_parton
     procedure :: is_gauge => field_data_is_gauge
     procedure :: is_left_handed => field_data_is_left_handed
     procedure :: is_right_handed => field_data_is_right_handed
     procedure :: has_antiparticle => field_data_has_antiparticle
     procedure :: is_stable => field_data_is_stable
     procedure :: get_decays => field_data_get_decays
     procedure :: decays_isotropically => field_data_decays_isotropically
     procedure :: decays_diagonal => field_data_decays_diagonal
     procedure :: has_decay_helicity => field_data_has_decay_helicity
     procedure :: decay_helicity => field_data_decay_helicity
     procedure :: is_polarized => field_data_is_polarized
     procedure :: get_longname => field_data_get_longname
     procedure :: get_name => field_data_get_name
     procedure :: get_name_array => field_data_get_name_array
     procedure :: get_tex_name => field_data_get_tex_name
     procedure, private :: matches_name => field_data_matches_name
     procedure :: get_spin_type => field_data_get_spin_type
     procedure :: get_multiplicity => field_data_get_multiplicity
     procedure :: get_isospin_type => field_data_get_isospin_type
     procedure :: get_charge_type => field_data_get_charge_type
     procedure :: get_color_type => field_data_get_color_type
     procedure :: get_charge => field_data_get_charge
     procedure :: get_isospin => field_data_get_isospin
     procedure :: get_mass => field_data_get_mass
     procedure :: get_mass_sign => field_data_get_mass_sign
     procedure :: get_width => field_data_get_width
  end type field_data_t

  type :: field_data_p
     type(field_data_t), pointer :: p => null ()
  end type field_data_p

  type :: vertex_t
     private
     logical :: trilinear
     integer, dimension(:), allocatable :: pdg
     type(field_data_p), dimension(:), allocatable :: prt
   contains
     procedure :: write => vertex_write
     procedure :: init => vertex_init
     procedure :: copy_from => vertex_copy_from
     procedure :: get_match => vertex_get_match
  end type vertex_t

  type :: vertex_iterator_t
     private
     class(model_data_t), pointer :: model => null ()
     integer, dimension(:), allocatable :: pdg
     integer :: vertex_index = 0
     integer :: pdg_index = 0
     logical :: save_pdg_index
   contains
     procedure :: init => vertex_iterator_init
     procedure :: get_next_match => vertex_iterator_get_next_match
  end type vertex_iterator_t

  type :: vertex_table_entry_t
     private
     integer :: pdg1 = 0, pdg2 = 0
     integer :: n = 0
     integer, dimension(:), allocatable :: pdg3
  end type vertex_table_entry_t

  type :: vertex_table_t
     type(vertex_table_entry_t), dimension(:), allocatable :: entry
     integer :: n_collisions = 0
     integer(i32) :: mask
   contains
     procedure :: write => vertex_table_write
     procedure :: init => vertex_table_init
     procedure :: match => vertex_table_match
     procedure :: check => vertex_table_check
  end type vertex_table_t

  type :: model_data_t
     private
     type(string_t) :: name
     integer :: scheme = 0
     type(modelpar_real_t), dimension(:), pointer :: par_real => null ()
     type(modelpar_complex_t), dimension(:), pointer :: par_complex => null ()
     type(field_data_t), dimension(:), allocatable :: field
     type(vertex_t), dimension(:), allocatable :: vtx
     type(vertex_table_t) :: vt
   contains
     procedure :: final => model_data_final
     procedure :: write => model_data_write
     generic :: init => model_data_init
     procedure, private :: model_data_init
     procedure :: set_scheme_num => model_data_set_scheme_num
     procedure :: freeze_fields => model_data_freeze_fields
     procedure :: copy_from => model_data_copy
     procedure :: get_name => model_data_get_name
     procedure :: get_scheme_num => model_data_get_scheme_num
     procedure :: get_parameters_md5sum => model_data_get_parameters_md5sum
     procedure :: get_md5sum => model_data_get_md5sum
     generic :: init_par => model_data_init_par_real, model_data_init_par_complex
     procedure, private :: model_data_init_par_real
     procedure, private :: model_data_init_par_complex
     procedure :: get_n_real => model_data_get_n_real
     procedure :: get_n_complex => model_data_get_n_complex
     procedure :: real_parameters_to_array &
          => model_data_real_par_to_array
     procedure :: complex_parameters_to_array &
          => model_data_complex_par_to_array
     procedure :: real_parameters_from_array &
          => model_data_real_par_from_array
     procedure :: complex_parameters_from_array &
          => model_data_complex_par_from_array
     procedure :: real_parameters_to_c_array &
          => model_data_real_par_to_c_array
     procedure :: real_parameters_from_c_array &
          => model_data_real_par_from_c_array
     procedure :: get_par_real_ptr => model_data_get_par_real_ptr_index
     procedure :: get_par_complex_ptr => model_data_get_par_complex_ptr_index
     procedure :: get_par_data_ptr => model_data_get_par_data_ptr_name
     procedure :: get_real => model_data_get_par_real_value
     procedure :: get_complex => model_data_get_par_complex_value
     generic :: set_par => model_data_set_par_real, model_data_set_par_complex
     procedure, private :: model_data_set_par_real
     procedure, private :: model_data_set_par_complex
     procedure :: write_fields => model_data_write_fields
     procedure :: get_n_field => model_data_get_n_field
     generic :: get_pdg => &
          model_data_get_field_pdg_index, &
          model_data_get_field_pdg_name
     procedure, private :: model_data_get_field_pdg_index
     procedure, private :: model_data_get_field_pdg_name
     procedure :: get_all_pdg => model_data_get_all_pdg
     procedure :: get_field_array_ptr => model_data_get_field_array_ptr
     generic :: get_field_ptr => &
          model_data_get_field_ptr_name, &
          model_data_get_field_ptr_pdg
     procedure, private :: model_data_get_field_ptr_name
     procedure, private :: model_data_get_field_ptr_pdg
     procedure :: get_field_ptr_by_index => model_data_get_field_ptr_index
     procedure :: test_field => model_data_test_field_pdg
     procedure :: field_error => model_data_field_error
     procedure :: set_field_mass => model_data_set_field_mass_pdg
     procedure :: set_field_width => model_data_set_field_width_pdg
     procedure :: set_unstable => model_data_set_unstable
     procedure :: set_stable => model_data_set_stable
     procedure :: set_polarized => model_data_set_polarized
     procedure :: set_unpolarized => model_data_set_unpolarized
     procedure :: clear_unstable => model_clear_unstable
     procedure :: clear_polarized => model_clear_polarized
     procedure :: write_vertices => model_data_write_vertices
     generic :: set_vertex => &
          model_data_set_vertex_pdg, model_data_set_vertex_names
     procedure, private :: model_data_set_vertex_pdg
     procedure, private :: model_data_set_vertex_names
     procedure :: freeze_vertices => model_data_freeze_vertices
     procedure :: get_n_vtx => model_data_get_n_vtx
     procedure :: match_vertex => model_data_match_vertex
     procedure :: check_vertex => model_data_check_vertex
     procedure :: init_test => model_data_init_test
     procedure :: init_qed_test => model_data_init_qed_test
     procedure :: init_sm_test => model_data_init_sm_test
  end type model_data_t


  interface
    module subroutine par_write (par, unit)
      class(modelpar_data_t), intent(in) :: par
      integer, intent(in), optional :: unit
    end subroutine
    module subroutine par_show (par, l, u)
      class(modelpar_data_t), intent(in) :: par
      integer, intent(in) :: l, u
    end subroutine par_show
    module subroutine modelpar_data_init_real (par, name, value)
      class(modelpar_data_t), intent(out) :: par
      type(string_t), intent(in) :: name
      real(default), intent(in) :: value
    end subroutine modelpar_data_init_real
    module subroutine modelpar_data_init_complex (par, name, value)
      class(modelpar_data_t), intent(out) :: par
      type(string_t), intent(in) :: name
      complex(default), intent(in) :: value
    end subroutine modelpar_data_init_complex
    elemental module subroutine modelpar_data_set_real (par, value)
      class(modelpar_data_t), intent(inout) :: par
      real(default), intent(in) :: value
    end subroutine modelpar_data_set_real
    elemental module subroutine modelpar_data_set_complex (par, value)
      class(modelpar_data_t), intent(inout) :: par
      complex(default), intent(in) :: value
    end subroutine modelpar_data_set_complex
    module function modelpar_data_get_name (par) result (name)
      class(modelpar_data_t), intent(in) :: par
      type(string_t) :: name
    end function modelpar_data_get_name
    elemental module function modelpar_data_get_real (par) result (value)
      class(modelpar_data_t), intent(in), target :: par
      real(default) :: value
    end function modelpar_data_get_real
    elemental module function modelpar_data_get_complex (par) result (value)
      class(modelpar_data_t), intent(in), target :: par
      complex(default) :: value
    end function modelpar_data_get_complex
    module function modelpar_data_get_real_ptr (par) result (ptr)
      class(modelpar_data_t), intent(in), target :: par
      real(default), pointer :: ptr
    end function modelpar_data_get_real_ptr
    module function modelpar_data_get_complex_ptr (par) result (ptr)
      class(modelpar_data_t), intent(in), target :: par
      complex(default), pointer :: ptr
    end function modelpar_data_get_complex_ptr
    module subroutine field_data_init (prt, longname, pdg)
      class(field_data_t), intent(out) :: prt
      type(string_t), intent(in) :: longname
      integer, intent(in) :: pdg
    end subroutine field_data_init
    module subroutine field_data_copy_from (prt, prt_src)
      class(field_data_t), intent(inout) :: prt
      class(field_data_t), intent(in) :: prt_src
    end subroutine field_data_copy_from
    module subroutine field_data_set (prt, &
         is_visible, is_parton, is_gauge, is_left_handed, is_right_handed, &
         p_is_stable, p_decays_isotropically, p_decays_diagonal, &
         p_decay_helicity, &
         a_is_stable, a_decays_isotropically, a_decays_diagonal, &
         a_decay_helicity, &
         p_polarized, a_polarized, &
         name, anti, tex_name, tex_anti, &
         spin_type, isospin_type, charge_type, color_type, &
         mass_data, width_data, &
         p_decay, a_decay)
      class(field_data_t), intent(inout) :: prt
      logical, intent(in), optional :: is_visible, is_parton, is_gauge
      logical, intent(in), optional :: is_left_handed, is_right_handed
      logical, intent(in), optional :: p_is_stable
      logical, intent(in), optional :: p_decays_isotropically, p_decays_diagonal
      integer, intent(in), optional :: p_decay_helicity
      logical, intent(in), optional :: a_is_stable
      logical, intent(in), optional :: a_decays_isotropically, a_decays_diagonal
      integer, intent(in), optional :: a_decay_helicity
      logical, intent(in), optional :: p_polarized, a_polarized
      type(string_t), dimension(:), intent(in), optional :: name, anti
      type(string_t), intent(in), optional :: tex_name, tex_anti
      integer, intent(in), optional :: spin_type, isospin_type
      integer, intent(in), optional :: charge_type, color_type
      class(modelpar_data_t), intent(in), pointer, optional :: mass_data, width_data
      type(string_t), dimension(:), intent(in), optional :: p_decay, a_decay
    end subroutine field_data_set
    module subroutine field_data_set_multiplicity (prt)
      class(field_data_t), intent(inout) :: prt
    end subroutine field_data_set_multiplicity
    module subroutine field_data_set_mass (prt, mass)
      class(field_data_t), intent(inout) :: prt
      real(default), intent(in) :: mass
    end subroutine field_data_set_mass
    module subroutine field_data_set_width (prt, width)
      class(field_data_t), intent(inout) :: prt
      real(default), intent(in) :: width
    end subroutine field_data_set_width
    elemental module subroutine field_data_freeze (prt)
      class(field_data_t), intent(inout) :: prt
    end subroutine field_data_freeze
    module subroutine field_data_write (prt, unit)
      class(field_data_t), intent(in) :: prt
      integer, intent(in), optional :: unit
    end subroutine field_data_write
    module subroutine field_data_write_decays (prt, unit)
      class(field_data_t), intent(in) :: prt
      integer, intent(in), optional :: unit
    end subroutine field_data_write_decays
    module subroutine field_data_show (prt, l, u)
      class(field_data_t), intent(in) :: prt
      integer, intent(in) :: l, u
    end subroutine field_data_show
    elemental module function field_data_get_pdg (prt) result (pdg)
      integer :: pdg
      class(field_data_t), intent(in) :: prt
    end function field_data_get_pdg
    elemental module function field_data_get_pdg_anti (prt) result (pdg)
      integer :: pdg
      class(field_data_t), intent(in) :: prt
    end function field_data_get_pdg_anti
    elemental module function field_data_is_visible (prt) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
    end function field_data_is_visible
    elemental module function field_data_is_parton (prt) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
    end function field_data_is_parton
    elemental module function field_data_is_gauge (prt) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
    end function field_data_is_gauge
    elemental module function field_data_is_left_handed (prt) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
    end function field_data_is_left_handed
    elemental module function field_data_is_right_handed (prt) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
    end function field_data_is_right_handed
    elemental module function field_data_has_antiparticle (prt) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
    end function field_data_has_antiparticle
    elemental module function field_data_is_stable (prt, anti) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
      logical, intent(in), optional :: anti
    end function field_data_is_stable
    module subroutine field_data_get_decays (prt, decay, anti)
      class(field_data_t), intent(in) :: prt
      type(string_t), dimension(:), intent(out), allocatable :: decay
      logical, intent(in), optional :: anti
    end subroutine field_data_get_decays
    elemental module function field_data_decays_isotropically &
         (prt, anti) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
      logical, intent(in), optional :: anti
    end function field_data_decays_isotropically
    elemental module function field_data_decays_diagonal &
         (prt, anti) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
      logical, intent(in), optional :: anti
    end function field_data_decays_diagonal
    elemental module function field_data_has_decay_helicity &
         (prt, anti) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
      logical, intent(in), optional :: anti
    end function field_data_has_decay_helicity
    elemental module function field_data_decay_helicity &
         (prt, anti) result (hel)
      integer :: hel
      class(field_data_t), intent(in) :: prt
      logical, intent(in), optional :: anti
    end function field_data_decay_helicity
    elemental module function field_data_is_polarized (prt, anti) result (flag)
      logical :: flag
      class(field_data_t), intent(in) :: prt
      logical, intent(in), optional :: anti
    end function field_data_is_polarized
    pure module function field_data_get_longname (prt) result (name)
      type(string_t) :: name
      class(field_data_t), intent(in) :: prt
    end function field_data_get_longname
    pure module function field_data_get_name &
         (prt, is_antiparticle) result (name)
      type(string_t) :: name
      class(field_data_t), intent(in) :: prt
      logical, intent(in) :: is_antiparticle
    end function field_data_get_name
    module subroutine field_data_get_name_array (prt, is_antiparticle, name)
      class(field_data_t), intent(in) :: prt
      logical, intent(in) :: is_antiparticle
      type(string_t), dimension(:), allocatable, intent(inout) :: name
    end subroutine field_data_get_name_array
    elemental module function field_data_get_tex_name &
         (prt, is_antiparticle) result (name)
      type(string_t) :: name
      class(field_data_t), intent(in) :: prt
      logical, intent(in) :: is_antiparticle
    end function field_data_get_tex_name
    module function field_data_matches_name &
         (field, name, is_antiparticle) result (flag)
      class(field_data_t), intent(in) :: field
      type(string_t), intent(in) :: name
      logical, intent(in) :: is_antiparticle
      logical :: flag
    end function field_data_matches_name
    elemental module function field_data_get_spin_type (prt) result (type)
      integer :: type
      class(field_data_t), intent(in) :: prt
    end function field_data_get_spin_type
    elemental module function field_data_get_multiplicity (prt) result (type)
      integer :: type
      class(field_data_t), intent(in) :: prt
    end function field_data_get_multiplicity
    elemental module function field_data_get_isospin_type (prt) result (type)
      integer :: type
      class(field_data_t), intent(in) :: prt
    end function field_data_get_isospin_type
    elemental module function field_data_get_charge_type (prt) result (type)
      integer :: type
      class(field_data_t), intent(in) :: prt
    end function field_data_get_charge_type
    elemental module function field_data_get_color_type (prt) result (type)
      integer :: type
      class(field_data_t), intent(in) :: prt
    end function field_data_get_color_type
    elemental module function field_data_get_charge (prt) result (charge)
      real(default) :: charge
      class(field_data_t), intent(in) :: prt
    end function field_data_get_charge
    elemental module function field_data_get_isospin (prt) result (isospin)
      real(default) :: isospin
      class(field_data_t), intent(in) :: prt
    end function field_data_get_isospin
    elemental module function field_data_get_mass (prt) result (mass)
      real(default) :: mass
      class(field_data_t), intent(in) :: prt
    end function field_data_get_mass
    elemental module function field_data_get_mass_sign (prt) result (sgn)
      integer :: sgn
      class(field_data_t), intent(in) :: prt
    end function field_data_get_mass_sign
    elemental module function field_data_get_width (prt) result (width)
      real(default) :: width
      class(field_data_t), intent(in) :: prt
    end function field_data_get_width
    module subroutine find_model (model, PDG, model_A, model_B)
      class(model_data_t), pointer, intent(out) :: model
      integer, intent(in) :: PDG
      class(model_data_t), intent(in), target :: model_A, model_B
    end subroutine find_model
    module subroutine vertex_write (vtx, unit)
      class(vertex_t), intent(in) :: vtx
      integer, intent(in), optional :: unit
    end subroutine vertex_write
    module subroutine vertex_init (vtx, pdg, model)
      class(vertex_t), intent(out) :: vtx
      integer, dimension(:), intent(in) :: pdg
      type(model_data_t), intent(in), target, optional :: model
    end subroutine vertex_init
    module subroutine vertex_copy_from (vtx, old_vtx, new_model)
      class(vertex_t), intent(out) :: vtx
      class(vertex_t), intent(in) :: old_vtx
      type(model_data_t), intent(in), target, optional :: new_model
    end subroutine vertex_copy_from
    module subroutine vertex_get_match (vtx, pdg1, pdg2)
      class(vertex_t), intent(in) :: vtx
      integer, intent(in) :: pdg1
      integer, dimension(:), allocatable, intent(out) :: pdg2
    end subroutine vertex_get_match
    module subroutine vertex_iterator_init (it, model, pdg, save_pdg_index)
      class(vertex_iterator_t), intent(out) :: it
      class(model_data_t), intent(in), target :: model
      integer, dimension(:), intent(in) :: pdg
      logical, intent(in) :: save_pdg_index
    end subroutine vertex_iterator_init
    module subroutine vertex_iterator_get_next_match (it, pdg_match)
      class(vertex_iterator_t), intent(inout) :: it
      integer, dimension(:), allocatable, intent(out) :: pdg_match
    end subroutine vertex_iterator_get_next_match
    module subroutine vertex_table_write (vt, unit)
      class(vertex_table_t), intent(in) :: vt
      integer, intent(in), optional :: unit
    end subroutine vertex_table_write
    module subroutine vertex_table_init (vt, prt, vtx)
      class(vertex_table_t), intent(out) :: vt
      type(field_data_t), dimension(:), intent(in) :: prt
      type(vertex_t), dimension(:), intent(in) :: vtx
    end subroutine vertex_table_init
    module subroutine vertex_table_match (vt, pdg1, pdg2, pdg3)
      class(vertex_table_t), intent(in) :: vt
      integer, intent(in) :: pdg1, pdg2
      integer, dimension(:), allocatable, intent(out) :: pdg3
    end subroutine vertex_table_match
    module function vertex_table_check (vt, pdg1, pdg2, pdg3) result (flag)
      class(vertex_table_t), intent(in) :: vt
      integer, intent(in) :: pdg1, pdg2, pdg3
      logical :: flag
    end function vertex_table_check
    module subroutine model_data_final (model)
      class(model_data_t), intent(inout) :: model
    end subroutine model_data_final
    module subroutine model_data_write (model, unit, verbose, &
         show_md5sum, show_variables, show_parameters, &
         show_particles, show_vertices, show_scheme)
      class(model_data_t), intent(in) :: model
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
      logical, intent(in), optional :: show_md5sum
      logical, intent(in), optional :: show_variables
      logical, intent(in), optional :: show_parameters
      logical, intent(in), optional :: show_particles
      logical, intent(in), optional :: show_vertices
      logical, intent(in), optional :: show_scheme
    end subroutine model_data_write
    module subroutine model_data_init (model, name, &
         n_par_real, n_par_complex, n_field, n_vtx)
      class(model_data_t), intent(out) :: model
      type(string_t), intent(in) :: name
      integer, intent(in) :: n_par_real, n_par_complex
      integer, intent(in) :: n_field
      integer, intent(in) :: n_vtx
    end subroutine model_data_init
    module subroutine model_data_set_scheme_num (model, scheme)
      class(model_data_t), intent(inout) :: model
      integer, intent(in) :: scheme
    end subroutine model_data_set_scheme_num
    module subroutine model_data_freeze_fields (model)
      class(model_data_t), intent(inout) :: model
    end subroutine model_data_freeze_fields
    module subroutine model_data_copy (model, src)
      class(model_data_t), intent(inout), target :: model
      class(model_data_t), intent(in), target :: src
    end subroutine model_data_copy
    module function model_data_get_name (model) result (name)
      class(model_data_t), intent(in) :: model
      type(string_t) :: name
    end function model_data_get_name
    module function model_data_get_scheme_num (model) result (scheme)
      class(model_data_t), intent(in) :: model
      integer :: scheme
    end function model_data_get_scheme_num
    module function model_data_get_parameters_md5sum (model) result (par_md5sum)
      character(32) :: par_md5sum
      class(model_data_t), intent(in) :: model
    end function model_data_get_parameters_md5sum
    module function model_data_get_md5sum (model) result (md5sum)
      class(model_data_t), intent(in) :: model
      character(32) :: md5sum
    end function model_data_get_md5sum
    module subroutine model_data_init_par_real (model, i, name, value)
      class(model_data_t), intent(inout) :: model
      integer, intent(in) :: i
      type(string_t), intent(in) :: name
      real(default), intent(in) :: value
    end subroutine model_data_init_par_real
    module subroutine model_data_init_par_complex (model, i, name, value)
      class(model_data_t), intent(inout) :: model
      integer, intent(in) :: i
      type(string_t), intent(in) :: name
      complex(default), intent(in) :: value
    end subroutine model_data_init_par_complex
    module function model_data_get_n_real (model) result (n)
      class(model_data_t), intent(in) :: model
      integer :: n
    end function model_data_get_n_real
    module function model_data_get_n_complex (model) result (n)
      class(model_data_t), intent(in) :: model
      integer :: n
    end function model_data_get_n_complex
    module subroutine model_data_real_par_to_array (model, array)
      class(model_data_t), intent(in) :: model
      real(default), dimension(:), intent(inout) :: array
    end subroutine model_data_real_par_to_array
    module subroutine model_data_complex_par_to_array (model, array)
      class(model_data_t), intent(in) :: model
      complex(default), dimension(:), intent(inout) :: array
    end subroutine model_data_complex_par_to_array
    module subroutine model_data_real_par_from_array (model, array)
      class(model_data_t), intent(inout) :: model
      real(default), dimension(:), intent(in) :: array
    end subroutine model_data_real_par_from_array
    module subroutine model_data_complex_par_from_array (model, array)
      class(model_data_t), intent(inout) :: model
      complex(default), dimension(:), intent(in) :: array
    end subroutine model_data_complex_par_from_array
    module subroutine model_data_real_par_to_c_array (model, array)
      class(model_data_t), intent(in) :: model
      real(c_default_float), dimension(:), intent(inout) :: array
    end subroutine model_data_real_par_to_c_array
    module subroutine model_data_real_par_from_c_array (model, array)
      class(model_data_t), intent(inout) :: model
      real(c_default_float), dimension(:), intent(in) :: array
    end subroutine model_data_real_par_from_c_array
    module function model_data_get_par_real_ptr_index (model, i) result (ptr)
      class(model_data_t), intent(inout) :: model
      integer, intent(in) :: i
      class(modelpar_data_t), pointer :: ptr
    end function model_data_get_par_real_ptr_index
    module function model_data_get_par_complex_ptr_index (model, i) result (ptr)
      class(model_data_t), intent(inout) :: model
      integer, intent(in) :: i
      class(modelpar_data_t), pointer :: ptr
    end function model_data_get_par_complex_ptr_index
    module function model_data_get_par_data_ptr_name (model, name) result (ptr)
      class(model_data_t), intent(in) :: model
      type(string_t), intent(in) :: name
      class(modelpar_data_t), pointer :: ptr
    end function model_data_get_par_data_ptr_name
    module function model_data_get_par_real_value (model, name) result (value)
      class(model_data_t), intent(in) :: model
      type(string_t), intent(in) :: name
      real(default) :: value
    end function model_data_get_par_real_value
    module function model_data_get_par_complex_value &
         (model, name) result (value)
      class(model_data_t), intent(in) :: model
      type(string_t), intent(in) :: name
      complex(default) :: value
    end function model_data_get_par_complex_value
    module subroutine model_data_set_par_real (model, name, value)
      class(model_data_t), intent(inout) :: model
      type(string_t), intent(in) :: name
      real(default), intent(in) :: value
    end subroutine model_data_set_par_real
    module subroutine model_data_set_par_complex (model, name, value)
      class(model_data_t), intent(inout) :: model
      type(string_t), intent(in) :: name
      complex(default), intent(in) :: value
    end subroutine model_data_set_par_complex
    module subroutine model_data_write_fields (model, unit)
      class(model_data_t), intent(in) :: model
      integer, intent(in), optional :: unit
    end subroutine model_data_write_fields
    module function model_data_get_n_field (model) result (n)
      class(model_data_t), intent(in) :: model
      integer :: n
    end function model_data_get_n_field
    module function model_data_get_field_pdg_index (model, i) result (pdg)
      class(model_data_t), intent(in) :: model
      integer, intent(in) :: i
      integer :: pdg
    end function model_data_get_field_pdg_index
    module function model_data_get_field_pdg_name &
         (model, name, check) result (pdg)
      class(model_data_t), intent(in) :: model
      type(string_t), intent(in) :: name
      logical, intent(in), optional :: check
      integer :: pdg
    end function model_data_get_field_pdg_name
    module subroutine model_data_get_all_pdg (model, pdg)
      class(model_data_t), intent(in) :: model
      integer, dimension(:), allocatable, intent(inout) :: pdg
    end subroutine model_data_get_all_pdg
    module function model_data_get_field_array_ptr (model) result (ptr)
      class(model_data_t), intent(in), target :: model
      type(field_data_t), dimension(:), pointer :: ptr
    end function model_data_get_field_array_ptr
    module function model_data_get_field_ptr_name &
         (model, name, check) result (ptr)
      class(model_data_t), intent(in), target :: model
      type(string_t), intent(in) :: name
      logical, intent(in), optional :: check
      type(field_data_t), pointer :: ptr
    end function model_data_get_field_ptr_name
    module function model_data_get_field_ptr_pdg &
         (model, pdg, check) result (ptr)
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: pdg
      logical, intent(in), optional :: check
      type(field_data_t), pointer :: ptr
    end function model_data_get_field_ptr_pdg
    module function model_data_get_field_ptr_index (model, i) result (ptr)
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: i
      type(field_data_t), pointer :: ptr
    end function model_data_get_field_ptr_index
    module function model_data_test_field_pdg (model, pdg, check) result (exist)
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: pdg
      logical, intent(in), optional :: check
      logical :: exist
    end function model_data_test_field_pdg
    module subroutine model_data_field_error (model, check, name, pdg)
      class(model_data_t), intent(in) :: model
      logical, intent(in), optional :: check
      type(string_t), intent(in), optional :: name
      integer, intent(in), optional :: pdg
    end subroutine model_data_field_error
    module subroutine model_data_set_field_mass_pdg (model, pdg, value)
      class(model_data_t), intent(inout) :: model
      integer, intent(in) :: pdg
      real(default), intent(in) :: value
    end subroutine model_data_set_field_mass_pdg
    module subroutine model_data_set_field_width_pdg (model, pdg, value)
      class(model_data_t), intent(inout) :: model
      integer, intent(in) :: pdg
      real(default), intent(in) :: value
    end subroutine model_data_set_field_width_pdg
    module subroutine model_data_set_unstable &
         (model, pdg, decay, isotropic, diagonal, decay_helicity)
      class(model_data_t), intent(inout), target :: model
      integer, intent(in) :: pdg
      type(string_t), dimension(:), intent(in) :: decay
      logical, intent(in), optional :: isotropic, diagonal
      integer, intent(in), optional :: decay_helicity
    end subroutine model_data_set_unstable
    module subroutine model_data_set_stable (model, pdg)
      class(model_data_t), intent(inout), target :: model
      integer, intent(in) :: pdg
    end subroutine model_data_set_stable
    module subroutine model_data_set_polarized (model, pdg)
      class(model_data_t), intent(inout), target :: model
      integer, intent(in) :: pdg
    end subroutine model_data_set_polarized
    module subroutine model_data_set_unpolarized (model, pdg)
      class(model_data_t), intent(inout), target :: model
      integer, intent(in) :: pdg
    end subroutine model_data_set_unpolarized
    module subroutine model_clear_unstable (model)
      class(model_data_t), intent(inout), target :: model
    end subroutine model_clear_unstable
    module subroutine model_clear_polarized (model)
      class(model_data_t), intent(inout), target :: model
    end subroutine model_clear_polarized
    module subroutine model_data_write_vertices (model, unit, verbose)
      class(model_data_t), intent(in) :: model
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine model_data_write_vertices
    module subroutine model_data_set_vertex_pdg (model, i, pdg)
      class(model_data_t), intent(inout), target :: model
      integer, intent(in) :: i
      integer, dimension(:), intent(in) :: pdg
    end subroutine model_data_set_vertex_pdg
    module subroutine model_data_set_vertex_names (model, i, name)
      class(model_data_t), intent(inout), target :: model
      integer, intent(in) :: i
      type(string_t), dimension(:), intent(in) :: name
    end subroutine model_data_set_vertex_names
    module subroutine model_data_freeze_vertices (model)
      class(model_data_t), intent(inout) :: model
    end subroutine model_data_freeze_vertices
    module function model_data_get_n_vtx (model) result (n)
      class(model_data_t), intent(in) :: model
      integer :: n
    end function model_data_get_n_vtx
    module subroutine model_data_match_vertex (model, pdg1, pdg2, pdg3)
      class(model_data_t), intent(in) :: model
      integer, intent(in) :: pdg1, pdg2
      integer, dimension(:), allocatable, intent(out) :: pdg3
    end subroutine model_data_match_vertex
    module function model_data_check_vertex &
         (model, pdg1, pdg2, pdg3) result (flag)
      logical :: flag
      class(model_data_t), intent(in) :: model
      integer, intent(in) :: pdg1, pdg2, pdg3
    end function model_data_check_vertex
    module subroutine model_data_init_test (model)
      class(model_data_t), intent(out) :: model
    end subroutine model_data_init_test
    module subroutine model_data_init_qed_test (model)
      class(model_data_t), intent(out) :: model
    end subroutine model_data_init_qed_test
    module subroutine model_data_init_sm_test (model)
      class(model_data_t), intent(out) :: model
    end subroutine model_data_init_sm_test
  end interface

end module model_data
