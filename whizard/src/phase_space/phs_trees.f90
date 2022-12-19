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

module phs_trees

  use kinds, only: default
  use kinds, only: TC
  use iso_varying_string, string_t => varying_string
  use lorentz
  use permutations, only: permutation_t, permutation_size
  use permutations, only: permutation_init, permutation_find
  use permutations, only: tc_decay_level, tc_permute
  use model_data
  use flavors
  use resonances, only: resonance_history_t, resonance_info_t
  use mappings

  implicit none
  private

  public :: phs_prt_t
  public :: phs_tree_t
  public :: phs_tree_from_array
  public :: phs_tree_equivalent
  public :: phs_tree_find_msq_permutation
  public :: phs_tree_find_angle_permutation
  public :: phs_tree_combine_particles
  public :: phs_tree_setup_prt_combinations
  public :: phs_tree_set_momentum_links

  type :: phs_prt_t
     private
     logical :: defined = .false.
     type(vector4_t) :: p
     real(default) :: p2
   contains
     procedure :: set_defined => phs_prt_set_defined
     procedure :: set_undefined => phs_prt_set_undefined
     procedure :: set_momentum => phs_prt_set_momentum
     procedure :: set_msq => phs_prt_set_msq
     procedure :: is_defined => phs_prt_is_defined
     procedure :: get_momentum => phs_prt_get_momentum
     procedure :: get_msq => phs_prt_get_msq
     procedure :: combine => phs_prt_combine
     procedure :: write => phs_prt_write
     procedure :: check => phs_prt_check
  end type phs_prt_t

  type :: phs_branch_t
     private
     logical :: set = .false.
     logical :: inverted_decay = .false.
     logical :: inverted_axis = .false.
     integer(TC) :: mother = 0
     integer(TC) :: sibling = 0
     integer(TC) :: friend = 0
     integer(TC) :: origin = 0
     integer(TC), dimension(2) :: daughter = 0
     integer :: firstborn = 0
     logical :: has_children = .false.
     logical :: has_friend = .false.
     logical :: is_real = .false.
  end type phs_branch_t

  type :: phs_tree_t
     private
     integer :: n_branches, n_externals, n_in, n_msq, n_angles
     integer(TC) :: n_branches_tot, n_branches_out
     integer(TC) :: mask, mask_in, mask_out
     type(phs_branch_t), dimension(:), allocatable :: branch
     type(mapping_t), dimension(:), allocatable :: mapping
     real(default), dimension(:), allocatable :: mass_sum
     real(default), dimension(:), allocatable :: effective_mass
     real(default), dimension(:), allocatable :: effective_width
     logical :: real_phsp = .false.
     integer, dimension(:), allocatable :: momentum_link
   contains
     procedure :: init => phs_tree_init
     procedure :: final => phs_tree_final
     procedure :: write => phs_tree_write
     procedure :: from_array => phs_tree_from_array
     procedure :: flip_t_to_s_channel => phs_tree_flip_t_to_s_channel
     procedure :: canonicalize => phs_tree_canonicalize
     procedure :: init_mapping => phs_tree_init_mapping
     procedure :: set_mapping_parameters => phs_tree_set_mapping_parameters
     procedure :: assign_s_mapping => phs_tree_assign_s_mapping
     procedure :: set_mass_sum => phs_tree_set_mass_sum
     procedure :: set_effective_masses => phs_tree_set_effective_masses
     procedure :: set_step_mappings => phs_tree_set_step_mappings
     procedure :: extract_resonance_history => phs_tree_extract_resonance_history
     procedure :: compute_volume => phs_tree_compute_volume
     procedure :: compute_momenta_from_x => phs_tree_compute_momenta_from_x
     procedure :: compute_x_from_momenta => phs_tree_compute_x_from_momenta
     procedure :: reshuffle_mappings => phs_tree_reshuffle_mappings
  end type phs_tree_t


  interface
    elemental module subroutine phs_prt_set_defined (prt)
      class(phs_prt_t), intent(inout) :: prt
    end subroutine phs_prt_set_defined
    elemental module subroutine phs_prt_set_undefined (prt)
      class(phs_prt_t), intent(inout) :: prt
    end subroutine phs_prt_set_undefined
    elemental module subroutine phs_prt_set_momentum (prt, p)
      class(phs_prt_t), intent(inout) :: prt
      type(vector4_t), intent(in) :: p
    end subroutine phs_prt_set_momentum
    elemental module subroutine phs_prt_set_msq (prt, p2)
      class(phs_prt_t), intent(inout) :: prt
      real(default), intent(in) :: p2
    end subroutine phs_prt_set_msq
    elemental module function phs_prt_is_defined (prt) result (defined)
      logical :: defined
      class(phs_prt_t), intent(in) :: prt
    end function phs_prt_is_defined
    elemental module function phs_prt_get_momentum (prt) result (p)
      type(vector4_t) :: p
      class(phs_prt_t), intent(in) :: prt
    end function phs_prt_get_momentum
    elemental module function phs_prt_get_msq (prt) result (p2)
      real(default) :: p2
      class(phs_prt_t), intent(in) :: prt
    end function phs_prt_get_msq
    elemental module subroutine phs_prt_combine (prt, prt1, prt2)
      class(phs_prt_t), intent(inout) :: prt
      type(phs_prt_t), intent(in) :: prt1, prt2
    end subroutine phs_prt_combine
    module subroutine phs_prt_write (prt, unit)
      class(phs_prt_t), intent(in) :: prt
      integer, intent(in), optional :: unit
    end subroutine phs_prt_write
    elemental module subroutine phs_prt_check (prt)
      class(phs_prt_t), intent(inout) :: prt
    end subroutine phs_prt_check
    elemental module subroutine phs_tree_init &
         (tree, n_in, n_out, n_masses, n_angles)
      class(phs_tree_t), intent(inout) :: tree
      integer, intent(in) :: n_in, n_out, n_masses, n_angles
    end subroutine phs_tree_init
    elemental module subroutine phs_tree_final (tree)
      class(phs_tree_t), intent(inout) :: tree
    end subroutine phs_tree_final
    module subroutine phs_tree_write (tree, unit)
      class(phs_tree_t), intent(in) :: tree
      integer, intent(in), optional :: unit
    end subroutine phs_tree_write
    module subroutine phs_tree_from_array (tree, a)
      class(phs_tree_t), intent(inout) :: tree
      integer(TC), dimension(:), intent(in) :: a
    end subroutine phs_tree_from_array
    module subroutine phs_tree_flip_t_to_s_channel (tree)
      class(phs_tree_t), intent(inout) :: tree
    end subroutine phs_tree_flip_t_to_s_channel
    module subroutine phs_tree_canonicalize (tree)
      class(phs_tree_t), intent(inout) :: tree
    end subroutine phs_tree_canonicalize
    module subroutine phs_tree_init_mapping (tree, k, type, pdg, model)
      class(phs_tree_t), intent(inout) :: tree
      integer(TC), intent(in) :: k
      type(string_t), intent(in) :: type
      integer, intent(in) :: pdg
      class(model_data_t), intent(in), target :: model
    end subroutine phs_tree_init_mapping
    module subroutine phs_tree_set_mapping_parameters &
         (tree, mapping_defaults, variable_limits)
      class(phs_tree_t), intent(inout) :: tree
      type(mapping_defaults_t), intent(in) :: mapping_defaults
      logical, intent(in) :: variable_limits
    end subroutine phs_tree_set_mapping_parameters
    module subroutine phs_tree_assign_s_mapping (tree, mapping)
      class(phs_tree_t), intent(in) :: tree
      type(mapping_t), intent(out) :: mapping
    end subroutine phs_tree_assign_s_mapping
    module subroutine phs_tree_set_mass_sum (tree, flv)
      class(phs_tree_t), intent(inout) :: tree
      type(flavor_t), dimension(:), intent(in) :: flv
    end subroutine phs_tree_set_mass_sum
    module subroutine phs_tree_set_effective_masses (tree)
      class(phs_tree_t), intent(inout) :: tree
    end subroutine phs_tree_set_effective_masses
    module subroutine phs_tree_set_step_mappings &
         (tree, exp_type, variable_limits)
      class(phs_tree_t), intent(inout) :: tree
      logical, intent(in) :: exp_type
      logical, intent(in) :: variable_limits
    end subroutine phs_tree_set_step_mappings
    module subroutine phs_tree_extract_resonance_history (tree, res_history)
      class(phs_tree_t), intent(in) :: tree
      type(resonance_history_t), intent(out) :: res_history
    end subroutine phs_tree_extract_resonance_history
    module function phs_tree_equivalent (t1, t2, perm) result (is_equal)
      type(phs_tree_t), intent(in) :: t1, t2
      type(permutation_t), intent(in) :: perm
      logical :: equal, is_equal
    end function phs_tree_equivalent
    module subroutine phs_tree_find_msq_permutation &
         (tree1, tree2, perm2, msq_perm)
      type(phs_tree_t), intent(in) :: tree1, tree2
      type(permutation_t), intent(in) :: perm2
      type(permutation_t), intent(out) :: msq_perm
    end subroutine phs_tree_find_msq_permutation
    module subroutine phs_tree_find_angle_permutation &
         (tree1, tree2, perm2, angle_perm, sig2)
      type(phs_tree_t), intent(in) :: tree1, tree2
      type(permutation_t), intent(in) :: perm2
      type(permutation_t), intent(out) :: angle_perm
      logical, dimension(:), allocatable, intent(out) :: sig2
    end subroutine phs_tree_find_angle_permutation
    module subroutine phs_tree_compute_volume (tree, sqrts, volume)
      class(phs_tree_t), intent(in) :: tree
      real(default), intent(in) :: sqrts
      real(default), intent(out) :: volume
    end subroutine phs_tree_compute_volume
    module subroutine phs_tree_compute_momenta_from_x &
         (tree, prt, factor, volume, sqrts, x, ok)
      class(phs_tree_t), intent(inout) :: tree
      type(phs_prt_t), dimension(:), intent(inout) :: prt
      real(default), intent(out) :: factor, volume
      real(default), intent(in) :: sqrts
      real(default), dimension(:), intent(in) :: x
      logical, intent(out) :: ok
    end subroutine phs_tree_compute_momenta_from_x
    module subroutine phs_tree_compute_x_from_momenta &
         (tree, prt, factor, sqrts, x)
      class(phs_tree_t), intent(inout) :: tree
      type(phs_prt_t), dimension(:), intent(in) :: prt
      real(default), intent(out) :: factor
      real(default), intent(in) :: sqrts
      real(default), dimension(:), intent(inout) :: x
    end subroutine phs_tree_compute_x_from_momenta
    module subroutine phs_tree_combine_particles (tree, prt)
      type(phs_tree_t), intent(in) :: tree
      type(phs_prt_t), dimension(:), intent(inout) :: prt
    end subroutine phs_tree_combine_particles
    module subroutine phs_tree_setup_prt_combinations (tree, comb)
      type(phs_tree_t), intent(in) :: tree
      integer, dimension(:,:), intent(out) :: comb
    end subroutine phs_tree_setup_prt_combinations
    module subroutine phs_tree_reshuffle_mappings (tree)
     class(phs_tree_t), intent(inout) :: tree
    end subroutine phs_tree_reshuffle_mappings
    module subroutine phs_tree_set_momentum_links (tree, list)
      type(phs_tree_t), intent(inout) :: tree
      integer, dimension(:), allocatable :: list
    end subroutine phs_tree_set_momentum_links
  end interface

end module phs_trees
