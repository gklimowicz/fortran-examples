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

module pdg_arrays

  implicit none
  private

  public :: pdg_array_t
  public :: pdg_array_write_set
  public :: assignment(=)
  public :: operator(//)
  public :: operator(.match.)
  public :: is_quark
  public :: is_gluon
  public :: is_photon
  public :: is_colored
  public :: is_lepton
  public :: is_charged_lepton
  public :: is_fermion
  public :: is_massless_vector
  public :: is_massive_vector
  public :: is_vector
  public :: is_elementary
  public :: is_ew_boson_scalar
  public :: query_coupling_powers
  public :: qcd_ew_interferences
  public :: flv_eqv_expr_class
  public :: operator(<)
  public :: operator(>)
  public :: operator(<=)
  public :: operator(>=)
  public :: operator(==)
  public :: operator(/=)
  public :: operator(.eqv.)
  public :: operator(.neqv.)
  public :: sort_abs
  public :: pdg_list_t

  type :: pdg_array_t
     private
     integer, dimension(:), allocatable :: pdg
   contains
     procedure :: write => pdg_array_write
     procedure :: init => pdg_array_init
     procedure :: delete => pdg_array_delete
     procedure :: merge => pdg_array_merge
     procedure :: get_length => pdg_array_get_length
     procedure :: get => pdg_array_get
     procedure :: set => pdg_array_set
     procedure :: add => pdg_array_add
     procedure :: replace => pdg_array_replace
     procedure :: has_colored_particles => pdg_array_has_colored_particles
     procedure :: sort_abs => pdg_array_sort_abs
     procedure :: intersect => pdg_array_intersect
     procedure :: search_for_particle => pdg_array_search_for_particle
     procedure :: invert => pdg_array_invert
  end type pdg_array_t

  type :: pdg_list_t
     type(pdg_array_t), dimension(:), allocatable :: a
   contains
     procedure :: write => pdg_list_write
     generic :: init => pdg_list_init_size
     procedure, private :: pdg_list_init_size
     generic :: init => pdg_list_init_int_array
     procedure, private :: pdg_list_init_int_array
     generic :: set => pdg_list_set_int
     generic :: set => pdg_list_set_int_array
     generic :: set => pdg_list_set_pdg_array
     procedure, private :: pdg_list_set_int
     procedure, private :: pdg_list_set_int_array
     procedure, private :: pdg_list_set_pdg_array
     procedure :: get_size => pdg_list_get_size
     procedure :: get => pdg_list_get
     procedure :: is_regular => pdg_list_is_regular
     procedure :: sort_abs => pdg_list_sort_abs
     generic :: operator (==) => pdg_list_eq
     procedure, private :: pdg_list_eq
     generic :: operator (<) => pdg_list_lt
     procedure, private :: pdg_list_lt
     procedure :: replace => pdg_list_replace
     procedure :: fusion => pdg_list_fusion
     procedure :: get_pdg_sizes => pdg_list_get_pdg_sizes
     procedure :: match_replace => pdg_list_match_replace
     generic :: operator (.match.) => pdg_list_match_pdg_array
     procedure, private :: pdg_list_match_pdg_array
     procedure :: find_match => pdg_list_find_match_pdg_array
     procedure :: create_pdg_array => pdg_list_create_pdg_array
     procedure :: create_antiparticles => pdg_list_create_antiparticles
     procedure :: search_for_particle => pdg_list_search_for_particle
     procedure :: contains_colored_particles => pdg_list_contains_colored_particles
  end type pdg_list_t


  interface assignment(=)
     module procedure pdg_array_from_int_array
     module procedure pdg_array_from_int
     module procedure int_array_from_pdg_array
  end interface

  interface operator(//)
     module procedure concat_pdg_arrays
  end interface

  interface operator(.match.)
     module procedure pdg_array_match_integer
     module procedure pdg_array_match_pdg_array
  end interface

  interface operator(<)
     module procedure pdg_array_lt
  end interface
  interface operator(>)
     module procedure pdg_array_gt
  end interface
  interface operator(<=)
     module procedure pdg_array_le
  end interface
  interface operator(>=)
     module procedure pdg_array_ge
  end interface
  interface operator(==)
     module procedure pdg_array_eq
  end interface
  interface operator(/=)
     module procedure pdg_array_ne
  end interface

  interface operator(.eqv.)
     module procedure pdg_array_equivalent
  end interface
  interface operator(.neqv.)
     module procedure pdg_array_inequivalent
  end interface

  interface sort_abs
     module procedure pdg_array_sort_abs
  end interface


  interface
    module subroutine pdg_array_write (aval, unit)
      class(pdg_array_t), intent(in) :: aval
      integer, intent(in), optional :: unit
    end subroutine pdg_array_write
    module subroutine pdg_array_write_set (aval, unit)
      type(pdg_array_t), intent(in), dimension(:) :: aval
      integer, intent(in), optional :: unit
    end subroutine pdg_array_write_set
    module subroutine pdg_array_from_int_array (aval, iarray)
      type(pdg_array_t), intent(out) :: aval
      integer, dimension(:), intent(in) :: iarray
    end subroutine pdg_array_from_int_array
    elemental module subroutine pdg_array_from_int (aval, int)
      type(pdg_array_t), intent(out) :: aval
      integer, intent(in) :: int
    end subroutine pdg_array_from_int
    module subroutine int_array_from_pdg_array (iarray, aval)
      integer, dimension(:), allocatable, intent(out) :: iarray
      type(pdg_array_t), intent(in) :: aval
    end subroutine int_array_from_pdg_array
    module subroutine pdg_array_init (aval, n_elements)
      class(pdg_array_t), intent(inout) :: aval
      integer, intent(in) :: n_elements
    end subroutine pdg_array_init
    module subroutine pdg_array_delete (aval)
      class(pdg_array_t), intent(inout) :: aval
    end subroutine pdg_array_delete
    module subroutine pdg_array_merge (aval1, aval2)
      class(pdg_array_t), intent(inout) :: aval1
      type(pdg_array_t), intent(in) :: aval2
    end subroutine pdg_array_merge
    elemental module function pdg_array_get_length (aval) result (n)
      class(pdg_array_t), intent(in) :: aval
      integer :: n
    end function pdg_array_get_length
    elemental module function pdg_array_get (aval, i) result (pdg)
      class(pdg_array_t), intent(in) :: aval
      integer, intent(in), optional :: i
      integer :: pdg
    end function pdg_array_get
    module subroutine pdg_array_set (aval, i, pdg)
      class(pdg_array_t), intent(inout) :: aval
      integer, intent(in) :: i
      integer, intent(in) :: pdg
    end subroutine pdg_array_set
    module function pdg_array_add (aval, aval_add) result (aval_out)
      type(pdg_array_t) :: aval_out
      class(pdg_array_t), intent(in) :: aval
      type(pdg_array_t), intent(in) :: aval_add
    end function pdg_array_add
    module function pdg_array_replace (aval, i, pdg_new) result (aval_new)
      class(pdg_array_t), intent(in) :: aval
      integer, intent(in) :: i
      integer, dimension(:), intent(in) :: pdg_new
      type(pdg_array_t) :: aval_new
    end function pdg_array_replace
    module function concat_pdg_arrays (aval1, aval2) result (aval)
      type(pdg_array_t) :: aval
      type(pdg_array_t), intent(in) :: aval1, aval2
    end function concat_pdg_arrays
    elemental module function pdg_array_match_integer (aval, pdg) result (flag)
      logical :: flag
      type(pdg_array_t), intent(in) :: aval
      integer, intent(in) :: pdg
    end function pdg_array_match_integer
    elemental module function is_quark (pdg_nr)
      logical :: is_quark
      integer, intent(in) :: pdg_nr
    end function is_quark
    elemental module function is_gluon (pdg_nr)
      logical :: is_gluon
      integer, intent(in) :: pdg_nr
    end function is_gluon
    elemental module function is_photon (pdg_nr)
      logical :: is_photon
      integer, intent(in) :: pdg_nr
    end function is_photon
    elemental module function is_colored (pdg_nr)
      logical :: is_colored
      integer, intent(in) :: pdg_nr
    end function is_colored
    elemental module function is_lepton (pdg_nr)
      logical :: is_lepton
      integer, intent(in) :: pdg_nr
    end function is_lepton
    elemental module function is_charged_lepton (pdg_nr)
      logical :: is_charged_lepton
      integer, intent(in) :: pdg_nr
    end function is_charged_lepton
    elemental module function is_fermion (pdg_nr)
      logical :: is_fermion
      integer, intent(in) :: pdg_nr
    end function is_fermion
    elemental module function is_massless_vector (pdg_nr)
      integer, intent(in) :: pdg_nr
      logical :: is_massless_vector
    end function is_massless_vector
    elemental module function is_massive_vector (pdg_nr)
      integer, intent(in) :: pdg_nr
      logical :: is_massive_vector
    end function is_massive_vector
    elemental module function is_vector (pdg_nr)
      integer, intent(in) :: pdg_nr
      logical :: is_vector
    end function is_vector
    elemental module function is_elementary (pdg_nr)
      integer, intent(in) :: pdg_nr
      logical :: is_elementary
    end function is_elementary
    elemental module function is_ew_boson_scalar (pdg_nr)
      integer, intent(in) :: pdg_nr
      logical :: is_ew_boson_scalar
    end function is_ew_boson_scalar
    module function pdg_array_has_colored_particles (pdg) result (colored)
      class(pdg_array_t), intent(in) :: pdg
      logical :: colored
    end function pdg_array_has_colored_particles
    module function query_coupling_powers (flv, a_power, as_power) result (valid)
       integer, intent(in), dimension(:) :: flv
       integer, intent(in) :: a_power, as_power
       logical :: valid
    end function query_coupling_powers
    module function qcd_ew_interferences (flv) result (valid)
       integer, intent(in), dimension(:) :: flv
       logical :: valid
    end function qcd_ew_interferences
    module function flv_eqv_expr_class (flv) result (assign_qgA)
       integer, intent(in) :: flv
       logical, dimension(3) :: assign_qgA
    end function flv_eqv_expr_class
    module function pdg_array_match_pdg_array (aval1, aval2) result (flag)
      logical :: flag
      type(pdg_array_t), intent(in) :: aval1, aval2
    end function pdg_array_match_pdg_array
    elemental module function pdg_array_lt (aval1, aval2) result (flag)
      type(pdg_array_t), intent(in) :: aval1, aval2
      logical :: flag
    end function pdg_array_lt
    elemental module function pdg_array_gt (aval1, aval2) result (flag)
      type(pdg_array_t), intent(in) :: aval1, aval2
      logical :: flag
    end function pdg_array_gt
    elemental module function pdg_array_le (aval1, aval2) result (flag)
      type(pdg_array_t), intent(in) :: aval1, aval2
      logical :: flag
    end function pdg_array_le
    elemental module function pdg_array_ge (aval1, aval2) result (flag)
      type(pdg_array_t), intent(in) :: aval1, aval2
      logical :: flag
    end function pdg_array_ge
    elemental module function pdg_array_eq (aval1, aval2) result (flag)
      type(pdg_array_t), intent(in) :: aval1, aval2
      logical :: flag
    end function pdg_array_eq
    elemental module function pdg_array_ne (aval1, aval2) result (flag)
      type(pdg_array_t), intent(in) :: aval1, aval2
      logical :: flag
    end function pdg_array_ne
    elemental module function pdg_array_equivalent (aval1, aval2) result (eq)
      logical :: eq
      type(pdg_array_t), intent(in) :: aval1, aval2
    end function pdg_array_equivalent
    elemental module function pdg_array_inequivalent (aval1, aval2) result (neq)
      logical :: neq
      type(pdg_array_t), intent(in) :: aval1, aval2
    end function pdg_array_inequivalent
    module function pdg_array_sort_abs (aval1, unique) result (aval2)
      class(pdg_array_t), intent(in) :: aval1
      logical, intent(in), optional :: unique
      type(pdg_array_t) :: aval2
    end function pdg_array_sort_abs
    module function pdg_array_intersect (aval1, match) result (aval2)
     class(pdg_array_t), intent(in) :: aval1
     integer, dimension(:) :: match
     type(pdg_array_t) :: aval2
    end function pdg_array_intersect
    elemental module function pdg_array_search_for_particle (pdg, i_part) result (found)
      class(pdg_array_t), intent(in) :: pdg
      integer, intent(in) :: i_part
      logical :: found
    end function pdg_array_search_for_particle
    module function pdg_array_invert (pdg) result (pdg_inverse)
      class(pdg_array_t), intent(in) :: pdg
      type(pdg_array_t) :: pdg_inverse
    end function pdg_array_invert
    module subroutine pdg_list_write (object, unit)
      class(pdg_list_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine pdg_list_write
    module subroutine pdg_list_init_size (pl, n)
      class(pdg_list_t), intent(out) :: pl
      integer, intent(in) :: n
    end subroutine pdg_list_init_size
    module subroutine pdg_list_init_int_array (pl, pdg)
      class(pdg_list_t), intent(out) :: pl
      integer, dimension(:), intent(in) :: pdg
    end subroutine pdg_list_init_int_array
    module subroutine pdg_list_set_int (pl, i, pdg)
      class(pdg_list_t), intent(inout) :: pl
      integer, intent(in) :: i
      integer, intent(in) :: pdg
    end subroutine pdg_list_set_int
    module subroutine pdg_list_set_int_array (pl, i, pdg)
      class(pdg_list_t), intent(inout) :: pl
      integer, intent(in) :: i
      integer, dimension(:), intent(in) :: pdg
    end subroutine pdg_list_set_int_array
    module subroutine pdg_list_set_pdg_array (pl, i, pa)
      class(pdg_list_t), intent(inout) :: pl
      integer, intent(in) :: i
      type(pdg_array_t), intent(in) :: pa
    end subroutine pdg_list_set_pdg_array
    module function pdg_list_get_size (pl) result (n)
      class(pdg_list_t), intent(in) :: pl
      integer :: n
    end function pdg_list_get_size
    module function pdg_list_get (pl, i) result (pa)
      type(pdg_array_t) :: pa
      class(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: i
    end function pdg_list_get
    module function pdg_list_is_regular (pl) result (flag)
      class(pdg_list_t), intent(in) :: pl
      logical :: flag
    end function pdg_list_is_regular
    module function pdg_list_sort_abs (pl, n_in) result (pl_sorted)
      class(pdg_list_t), intent(in) :: pl
      integer, intent(in), optional :: n_in
      type(pdg_list_t) :: pl_sorted
    end function pdg_list_sort_abs
    module function pdg_list_eq (pl1, pl2) result (flag)
      class(pdg_list_t), intent(in) :: pl1, pl2
      logical :: flag
    end function pdg_list_eq
    module function pdg_list_lt (pl1, pl2) result (flag)
      class(pdg_list_t), intent(in) :: pl1, pl2
      logical :: flag
    end function pdg_list_lt
    module function pdg_list_replace (pl, i, pl_insert, n_in) result (pl_out)
      type(pdg_list_t) :: pl_out
      class(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: i
      class(pdg_list_t), intent(in) :: pl_insert
      integer, intent(in), optional :: n_in
    end function pdg_list_replace
    module function pdg_list_fusion (pl, pl_insert, i, check_if_existing) result (pl_out)
      type(pdg_list_t) :: pl_out
      class(pdg_list_t), intent(in) :: pl
      type(pdg_list_t), intent(in) :: pl_insert
      integer, intent(in) :: i
      logical, intent(in) :: check_if_existing
    end function pdg_list_fusion
    module function pdg_list_get_pdg_sizes (pl) result (i_size)
      integer, dimension(:), allocatable :: i_size
      class(pdg_list_t), intent(in) :: pl
    end function pdg_list_get_pdg_sizes
    module subroutine pdg_list_match_replace (pl, pl_match, success)
      class(pdg_list_t), intent(inout) :: pl
      class(pdg_list_t), intent(in) :: pl_match
      logical, intent(out) :: success
    end subroutine pdg_list_match_replace
    module function pdg_list_match_pdg_array (pl, pa) result (flag)
      class(pdg_list_t), intent(in) :: pl
      type(pdg_array_t), intent(in) :: pa
      logical :: flag
    end function pdg_list_match_pdg_array
    module function pdg_list_find_match_pdg_array (pl, pa, mask) result (i)
      class(pdg_list_t), intent(in) :: pl
      type(pdg_array_t), intent(in) :: pa
      logical, dimension(:), intent(in), optional :: mask
      integer :: i
    end function pdg_list_find_match_pdg_array
    module subroutine pdg_list_create_pdg_array (pl, pdg)
      class(pdg_list_t), intent(in) :: pl
      type(pdg_array_t), dimension(:), intent(inout), allocatable :: pdg
    end subroutine pdg_list_create_pdg_array
    module subroutine pdg_list_create_antiparticles (pl, pl_anti, n_new_particles)
      class(pdg_list_t), intent(in) :: pl
      type(pdg_list_t), intent(out) :: pl_anti
      integer, intent(out) :: n_new_particles
    end subroutine pdg_list_create_antiparticles
    elemental module function pdg_list_search_for_particle (pl, i_part) result (found)
      logical :: found
      class(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: i_part
    end function pdg_list_search_for_particle
    module function pdg_list_contains_colored_particles (pl) result (colored)
      class(pdg_list_t), intent(in) :: pl
      logical :: colored
    end function pdg_list_contains_colored_particles
  end interface

end module pdg_arrays
