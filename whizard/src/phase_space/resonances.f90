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

module resonances

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use model_data, only: model_data_t
  use flavors, only: flavor_t

  implicit none
  private

  public :: resonance_contributors_t
  public :: resonance_info_t
  public :: resonance_history_t
  public :: resonance_tree_t
  public :: resonance_history_set_t

  integer, parameter :: n_max_resonances = 10
  integer, parameter :: resonance_history_set_initial_size = 16

  type :: resonance_contributors_t
     integer, dimension(:), allocatable :: c
   contains
     procedure, private :: resonance_contributors_equal
     generic :: operator(==) => resonance_contributors_equal
     procedure, private :: resonance_contributors_assign
     generic :: assignment(=) => resonance_contributors_assign
  end type resonance_contributors_t

  type :: resonance_info_t
     type(flavor_t) :: flavor
     type(resonance_contributors_t) :: contributors
  contains
     procedure :: copy => resonance_info_copy
     procedure :: write => resonance_info_write
     procedure, private :: resonance_info_init_pdg
     procedure, private :: resonance_info_init_flv
     generic :: init => resonance_info_init_pdg, resonance_info_init_flv
     procedure, private :: resonance_info_equal
     generic :: operator(==) => resonance_info_equal
     procedure :: mapping => resonance_info_mapping
     procedure, private :: get_n_contributors => resonance_info_get_n_contributors
     procedure, private :: contains => resonance_info_contains
     procedure :: evaluate_distance => resonance_info_evaluate_distance
     procedure :: evaluate_gaussian => resonance_info_evaluate_gaussian
     procedure :: is_on_shell => resonance_info_is_on_shell
     procedure :: as_omega_string => resonance_info_as_omega_string
  end type resonance_info_t

  type :: resonance_history_t
     type(resonance_info_t), dimension(:), allocatable :: resonances
     integer :: n_resonances = 0
  contains
     procedure :: clear => resonance_history_clear
     procedure :: copy => resonance_history_copy
     procedure :: write => resonance_history_write
     procedure, private :: resonance_history_assign
     generic :: assignment(=) => resonance_history_assign
     procedure, private :: resonance_history_equal
     generic :: operator(==) => resonance_history_equal
     procedure, private :: resonance_history_contains
     generic :: operator(.contains.) => resonance_history_contains
     procedure :: add_resonance => resonance_history_add_resonance
     procedure :: remove_resonance => resonance_history_remove_resonance
     procedure :: add_offset => resonance_history_add_offset
     procedure :: contains_leg => resonance_history_contains_leg
     procedure :: mapping => resonance_history_mapping
     procedure :: only_has_n_contributors => &
          resonance_history_only_has_n_contributors
     procedure :: has_flavor => resonance_history_has_flavor
     procedure :: evaluate_distances => resonance_history_evaluate_distances
     procedure :: evaluate_gaussian => resonance_history_evaluate_gaussian
     procedure :: is_on_shell => resonance_history_is_on_shell
     procedure :: as_omega_string => resonance_history_as_omega_string
     procedure :: to_tree => resonance_history_to_tree
  end type resonance_history_t

  type :: resonance_branch_t
     integer :: i = 0
     type(flavor_t) :: flv
     integer, dimension(:), allocatable :: r_child
     integer, dimension(:), allocatable :: o_child
  end type resonance_branch_t

  type :: resonance_tree_t
     private
     integer :: n = 0
     type(resonance_branch_t), dimension(:), allocatable :: branch
   contains
     procedure :: write => resonance_tree_write
     procedure :: get_n_resonances => resonance_tree_get_n_resonances
     procedure :: get_flv => resonance_tree_get_flv
     procedure :: get_children => resonance_tree_get_children
  end type resonance_tree_t


  type :: index_array_t
     integer, dimension(:), allocatable :: i
  end type index_array_t

  type :: resonance_history_set_t
     private
     logical :: complete = .false.
     integer :: n_filter = 0
     type(resonance_history_t), dimension(:), allocatable :: history
     type(index_array_t), dimension(:), allocatable :: contains_this
     type(resonance_tree_t), dimension(:), allocatable :: tree
     integer :: last = 0
   contains
     procedure :: write => resonance_history_set_write
     procedure :: init => resonance_history_set_init
     procedure :: enter => resonance_history_set_enter
     procedure :: freeze => resonance_history_set_freeze
     procedure :: determine_on_shell_histories &
          => resonance_history_set_determine_on_shell_histories
     procedure :: evaluate_gaussian => resonance_history_set_evaluate_gaussian
     procedure :: get_n_history => resonance_history_set_get_n_history
     procedure :: get_history => resonance_history_set_get_history
     procedure :: to_array => resonance_history_set_to_array
     procedure :: get_tree => resonance_history_set_get_tree
     procedure, private :: expand => resonance_history_set_expand
  end type resonance_history_set_t


  interface
    elemental module function resonance_contributors_equal &
         (c1, c2) result (equal)
      logical :: equal
      class(resonance_contributors_t), intent(in) :: c1, c2
    end function resonance_contributors_equal
    pure module subroutine resonance_contributors_assign &
         (contributors_out, contributors_in)
      class(resonance_contributors_t), intent(inout) :: contributors_out
      class(resonance_contributors_t), intent(in) :: contributors_in
    end subroutine resonance_contributors_assign
    module subroutine resonance_info_copy (resonance_in, resonance_out)
      class(resonance_info_t), intent(in) :: resonance_in
      type(resonance_info_t), intent(out) :: resonance_out
    end subroutine resonance_info_copy
    module subroutine resonance_info_write (resonance, unit, verbose)
      class(resonance_info_t), intent(in) :: resonance
      integer, optional, intent(in) :: unit
      logical, optional, intent(in) :: verbose
    end subroutine resonance_info_write
    module subroutine resonance_info_init_pdg &
         (resonance, mom_id, pdg, model, n_out)
      class(resonance_info_t), intent(out) :: resonance
      integer, intent(in) :: mom_id
      integer, intent(in) :: pdg, n_out
      class(model_data_t), intent(in), target :: model
    end subroutine resonance_info_init_pdg
    module subroutine resonance_info_init_flv (resonance, mom_id, flv, n_out)
      class(resonance_info_t), intent(out) :: resonance
      integer, intent(in) :: mom_id
      type(flavor_t), intent(in) :: flv
      integer, intent(in) :: n_out
    end subroutine resonance_info_init_flv
    elemental module function resonance_info_equal (r1, r2) result (equal)
      logical :: equal
      class(resonance_info_t), intent(in) :: r1, r2
    end function resonance_info_equal
    module function resonance_info_mapping (resonance, s) result (bw)
      real(default) :: bw
      class(resonance_info_t), intent(in) :: resonance
      real(default), intent(in) :: s
    end function resonance_info_mapping
    elemental module function resonance_info_get_n_contributors &
         (resonance) result (n)
      class(resonance_info_t), intent(in) :: resonance
      integer :: n
    end function resonance_info_get_n_contributors
    elemental module function resonance_info_contains &
         (resonance, c) result (flag)
      class(resonance_info_t), intent(in) :: resonance
      integer, intent(in) :: c
      logical :: flag
    end function resonance_info_contains
    module subroutine resonance_history_clear (res_hist)
      class(resonance_history_t), intent(out) :: res_hist
    end subroutine resonance_history_clear
    module subroutine resonance_history_copy (res_hist_in, res_hist_out)
      class(resonance_history_t), intent(in) :: res_hist_in
      type(resonance_history_t), intent(out) :: res_hist_out
    end subroutine resonance_history_copy
    module subroutine resonance_history_write (res_hist, unit, verbose, indent)
      class(resonance_history_t), intent(in) :: res_hist
      integer, optional, intent(in) :: unit
      logical, optional, intent(in) :: verbose
      integer, optional, intent(in) :: indent
    end subroutine resonance_history_write
    module subroutine resonance_history_assign (res_hist_out, res_hist_in)
      class(resonance_history_t), intent(out) :: res_hist_out
      class(resonance_history_t), intent(in) :: res_hist_in
    end subroutine resonance_history_assign
    elemental module function resonance_history_equal (rh1, rh2) result (equal)
      logical :: equal
      class(resonance_history_t), intent(in) :: rh1, rh2
    end function resonance_history_equal
    elemental module function resonance_history_contains &
         (rh1, rh2) result (flag)
      logical :: flag
      class(resonance_history_t), intent(in) :: rh1, rh2
    end function resonance_history_contains
    module subroutine resonance_history_add_resonance (res_hist, resonance)
      class(resonance_history_t), intent(inout) :: res_hist
      type(resonance_info_t), intent(in) :: resonance
    end subroutine resonance_history_add_resonance
    module subroutine resonance_history_remove_resonance (res_hist, i_res)
      class(resonance_history_t), intent(inout) :: res_hist
      integer, intent(in) :: i_res
    end subroutine resonance_history_remove_resonance
    module subroutine resonance_history_add_offset (res_hist, n)
      class(resonance_history_t), intent(inout) :: res_hist
      integer, intent(in) :: n
    end subroutine resonance_history_add_offset
    module function resonance_history_contains_leg &
         (res_hist, i_leg) result (val)
      logical :: val
      class(resonance_history_t), intent(in) :: res_hist
      integer, intent(in) :: i_leg
    end function resonance_history_contains_leg
    module function resonance_history_mapping &
         (res_hist, p, i_gluon) result (p_map)
      real(default) :: p_map
      class(resonance_history_t), intent(in) :: res_hist
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in), optional :: i_gluon
    end function resonance_history_mapping
    module function resonance_history_only_has_n_contributors &
         (res_hist, n) result (value)
      logical :: value
      class(resonance_history_t), intent(in) :: res_hist
      integer, intent(in) :: n
    end function resonance_history_only_has_n_contributors
    module function resonance_history_has_flavor &
         (res_hist, flv) result (has_flv)
      logical :: has_flv
      class(resonance_history_t), intent(in) :: res_hist
      type(flavor_t), intent(in) :: flv
    end function resonance_history_has_flavor
    module subroutine resonance_info_evaluate_distance (res_info, p, dist)
      class(resonance_info_t), intent(in) :: res_info
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(out) :: dist
    end subroutine resonance_info_evaluate_distance
    module subroutine resonance_history_evaluate_distances (res_hist, p, dist)
      class(resonance_history_t), intent(in) :: res_hist
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), dimension(:), intent(out) :: dist
    end subroutine resonance_history_evaluate_distances
    module function resonance_info_evaluate_gaussian &
         (res_info, p, gw) result (factor)
      class(resonance_info_t), intent(in) :: res_info
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: gw
      real(default) :: factor
    end function resonance_info_evaluate_gaussian
    module function resonance_history_evaluate_gaussian &
         (res_hist, p, gw) result (factor)
      class(resonance_history_t), intent(in) :: res_hist
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: gw
      real(default) :: factor
    end function resonance_history_evaluate_gaussian
    module function resonance_info_is_on_shell (res_info, p, on_shell_limit) &
         result (flag)
      class(resonance_info_t), intent(in) :: res_info
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: on_shell_limit
      logical :: flag
    end function resonance_info_is_on_shell
    module function resonance_history_is_on_shell &
         (res_hist, p, on_shell_limit) result (flag)
      class(resonance_history_t), intent(in) :: res_hist
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: on_shell_limit
      logical :: flag
    end function resonance_history_is_on_shell
    module function resonance_info_as_omega_string &
         (res_info, n_in) result (string)
      class(resonance_info_t), intent(in) :: res_info
      integer, intent(in) :: n_in
      type(string_t) :: string
    end function resonance_info_as_omega_string
    module function resonance_history_as_omega_string &
         (res_hist, n_in) result (string)
      class(resonance_history_t), intent(in) :: res_hist
      integer, intent(in) :: n_in
      type(string_t) :: string
    end function resonance_history_as_omega_string
    module subroutine resonance_tree_write (tree, unit, indent)
      class(resonance_tree_t), intent(in) :: tree
      integer, intent(in), optional :: unit, indent
    end subroutine resonance_tree_write
    module function resonance_tree_get_n_resonances (tree) result (n)
      class(resonance_tree_t), intent(in) :: tree
      integer :: n
    end function resonance_tree_get_n_resonances
    module function resonance_tree_get_flv (tree, i) result (flv)
      class(resonance_tree_t), intent(in) :: tree
      integer, intent(in) :: i
      type(flavor_t) :: flv
    end function resonance_tree_get_flv
    module function resonance_tree_get_children (tree, i, offset_r, offset_o) &
         result (child)
      class(resonance_tree_t), intent(in) :: tree
      integer, intent(in) :: i, offset_r, offset_o
      integer, dimension(:), allocatable :: child
    end function resonance_tree_get_children
    module subroutine resonance_history_to_tree (res_hist, tree)
      class(resonance_history_t), intent(in) :: res_hist
      type(resonance_tree_t), intent(out) :: tree
    end subroutine resonance_history_to_tree
    module subroutine resonance_history_set_write &
         (res_set, unit, indent, show_trees)
      class(resonance_history_set_t), intent(in) :: res_set
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: indent
      logical, intent(in), optional :: show_trees
    end subroutine resonance_history_set_write
    module subroutine resonance_history_set_init &
         (res_set, n_filter, initial_size)
      class(resonance_history_set_t), intent(out) :: res_set
      integer, intent(in), optional :: n_filter
      integer, intent(in), optional :: initial_size
    end subroutine resonance_history_set_init
    module subroutine resonance_history_set_enter &
         (res_set, res_history, trivial)
      class(resonance_history_set_t), intent(inout) :: res_set
      type(resonance_history_t), intent(in) :: res_history
      logical, intent(in), optional :: trivial
    end subroutine resonance_history_set_enter
    module subroutine resonance_history_set_freeze (res_set)
      class(resonance_history_set_t), intent(inout) :: res_set
    end subroutine resonance_history_set_freeze
    module subroutine resonance_history_set_determine_on_shell_histories &
         (res_set, p, on_shell_limit, index_array)
      class(resonance_history_set_t), intent(in) :: res_set
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: on_shell_limit
      integer, dimension(:), allocatable, intent(out) :: index_array
    end subroutine resonance_history_set_determine_on_shell_histories
    module function resonance_history_set_evaluate_gaussian &
         (res_set, p, gw, i) result (factor)
      class(resonance_history_set_t), intent(in) :: res_set
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: gw
      integer, intent(in) :: i
      real(default) :: factor
    end function resonance_history_set_evaluate_gaussian
    module function resonance_history_set_get_n_history (res_set) result (n)
      class(resonance_history_set_t), intent(in) :: res_set
      integer :: n
    end function resonance_history_set_get_n_history
    module function resonance_history_set_get_history &
         (res_set, i) result (res_history)
      class(resonance_history_set_t), intent(in) :: res_set
      integer, intent(in) :: i
      type(resonance_history_t) :: res_history
    end function resonance_history_set_get_history
    module subroutine resonance_history_set_to_array (res_set, res_history)
      class(resonance_history_set_t), intent(in) :: res_set
      type(resonance_history_t), dimension(:), allocatable, intent(out) :: &
           res_history
    end subroutine resonance_history_set_to_array
    module subroutine resonance_history_set_get_tree (res_set, i, res_tree)
      class(resonance_history_set_t), intent(in) :: res_set
      integer, intent(in) :: i
      type(resonance_tree_t), intent(out) :: res_tree
    end subroutine resonance_history_set_get_tree
    module subroutine resonance_history_set_expand (res_set)
      class(resonance_history_set_t), intent(inout) :: res_set
    end subroutine resonance_history_set_expand
  end interface

end module resonances
