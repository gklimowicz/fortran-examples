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

module subevents

  use, intrinsic :: iso_c_binding !NODEP!

  use kinds, only: default
  use numeric_utils, only: pacify
  use c_particles
  use lorentz
  use pdg_arrays
  use jets

  implicit none
  private

  public :: prt_t
  public :: prt_init_combine
  public :: prt_get_pdg
  public :: prt_get_momentum
  public :: prt_get_msq
  public :: prt_is_polarized
  public :: prt_get_helicity
  public :: prt_is_colorized
  public :: prt_is_clustered
  public :: prt_is_recombinable
  public :: prt_is_photon
  public :: prt_is_parton
  public :: prt_is_lepton
  public :: prt_is_b_jet
  public :: prt_is_c_jet
  public :: prt_get_n_col
  public :: prt_get_n_acl
  public :: prt_get_color_indices
  public :: c_prt
  public :: prt_write
  public :: operator(.match.)
  public :: are_disjoint
  public :: subevt_t
  public :: subevt_init
  public :: subevt_polarize
  public :: subevt_colorize
  public :: subevt_join
  public :: subevt_combine
  public :: subevt_collect
  public :: subevt_cluster
  public :: subevt_recombine
  public :: subevt_select
  public :: subevt_extract
  public :: subevt_sort
  public :: subevt_select_pdg_code
  public :: pacify

  integer, parameter, public :: PRT_UNDEFINED = 0
  integer, parameter, public :: PRT_BEAM = -9
  integer, parameter, public :: PRT_INCOMING = 1
  integer, parameter, public :: PRT_OUTGOING = 2
  integer, parameter, public :: PRT_COMPOSITE = 3
  integer, parameter, public :: PRT_VIRTUAL = 4
  integer, parameter, public :: PRT_RESONANT = 5
  integer, parameter, public :: PRT_BEAM_REMNANT = 9


  type :: prt_t
     private
     integer :: type = PRT_UNDEFINED
     integer :: pdg
     logical :: polarized = .false.
     logical :: colorized = .false.
     logical :: clustered = .false.
     logical :: is_b_jet = .false.
     logical :: is_c_jet = .false.
     integer :: h
     type(vector4_t) :: p
     real(default) :: p2
     integer, dimension(:), allocatable :: src
     integer, dimension(:), allocatable :: col
     integer, dimension(:), allocatable :: acl
  end type prt_t

  type :: subevt_t
     private
     integer :: n_active = 0
     type(prt_t), dimension(:), allocatable :: prt
   contains
     procedure :: reset => subevt_reset
     procedure :: write => subevt_write
     procedure :: set_beam => subevt_set_beam
     procedure :: set_composite => subevt_set_composite
     procedure :: set_incoming => subevt_set_incoming
     procedure :: set_outgoing => subevt_set_outgoing
     procedure :: set_pdg_beam => subevt_set_pdg_beam
     procedure :: set_pdg_incoming => subevt_set_pdg_incoming
     procedure :: set_pdg_outgoing => subevt_set_pdg_outgoing
     procedure :: set_p_beam => subevt_set_p_beam
     procedure :: set_p_incoming => subevt_set_p_incoming
     procedure :: set_p_outgoing => subevt_set_p_outgoing
     procedure :: set_p2_beam => subevt_set_p2_beam
     procedure :: set_p2_incoming => subevt_set_p2_incoming
     procedure :: set_p2_outgoing => subevt_set_p2_outgoing
     procedure :: is_nonempty => subevt_is_nonempty
     procedure :: get_length => subevt_get_length
     procedure :: get_prt => subevt_get_prt
     procedure :: get_sqrts_hat => subevt_get_sqrts_hat
     procedure :: get_n_in => subevt_get_n_in
     procedure :: get_n_out => subevt_get_n_out
  end type subevt_t


  interface c_prt
     module procedure c_prt_from_prt
  end interface

  interface operator(.match.)
     module procedure prt_match
  end interface

  interface assignment(=)
     module procedure subevt_assign
  end interface

  interface c_prt
     module procedure c_prt_from_subevt
     module procedure c_prt_array_from_subevt
  end interface

  interface subevt_sort
     module procedure subevt_sort_pdg
     module procedure subevt_sort_int
     module procedure subevt_sort_real
  end interface

  interface pacify
     module procedure pacify_prt
     module procedure pacify_subevt
  end interface pacify


  interface
    module subroutine prt_init_combine (prt, prt1, prt2)
      type(prt_t), intent(out) :: prt
      type(prt_t), intent(in) :: prt1, prt2
    end subroutine prt_init_combine
    elemental module function prt_get_pdg (prt) result (pdg)
      integer :: pdg
      type(prt_t), intent(in) :: prt
    end function prt_get_pdg
    elemental module function prt_get_momentum (prt) result (p)
      type(vector4_t) :: p
      type(prt_t), intent(in) :: prt
    end function prt_get_momentum
    elemental module function prt_get_msq (prt) result (msq)
      real(default) :: msq
      type(prt_t), intent(in) :: prt
    end function prt_get_msq
    elemental module function prt_is_polarized (prt) result (flag)
      logical :: flag
      type(prt_t), intent(in) :: prt
    end function prt_is_polarized
    elemental module function prt_get_helicity (prt) result (h)
      integer :: h
      type(prt_t), intent(in) :: prt
    end function prt_get_helicity
    elemental module function prt_is_colorized (prt) result (flag)
      logical :: flag
      type(prt_t), intent(in) :: prt
    end function prt_is_colorized
    elemental module function prt_is_clustered (prt) result (flag)
      logical :: flag
      type(prt_t), intent(in) :: prt
    end function prt_is_clustered
    elemental module function prt_is_recombinable (prt) result (flag)
      logical :: flag
      type(prt_t), intent(in) :: prt
    end function prt_is_recombinable
    elemental module function prt_is_photon (prt) result (flag)
      logical :: flag
      type(prt_t), intent(in) :: prt
    end function prt_is_photon
    elemental module function prt_is_parton (prt) result (flag)
      logical :: flag
      type(prt_t), intent(in) :: prt
    end function prt_is_parton
    elemental module function prt_is_lepton (prt) result (flag)
      logical :: flag
      type(prt_t), intent(in) :: prt
    end function prt_is_lepton
    elemental module function prt_is_b_jet (prt) result (flag)
      logical :: flag
      type(prt_t), intent(in) :: prt
    end function prt_is_b_jet
    elemental module function prt_is_c_jet (prt) result (flag)
      logical :: flag
      type(prt_t), intent(in) :: prt
    end function prt_is_c_jet
    elemental module function prt_get_n_col (prt) result (n)
      integer :: n
      type(prt_t), intent(in) :: prt
    end function prt_get_n_col
    elemental module function prt_get_n_acl (prt) result (n)
      integer :: n
      type(prt_t), intent(in) :: prt
    end function prt_get_n_acl
    module subroutine prt_get_color_indices (prt, col, acl)
      type(prt_t), intent(in) :: prt
      integer, dimension(:), allocatable, intent(out) :: col, acl
    end subroutine prt_get_color_indices
    elemental module function c_prt_from_prt (prt) result (c_prt)
      type(c_prt_t) :: c_prt
      type(prt_t), intent(in) :: prt
    end function c_prt_from_prt
    module subroutine prt_write (prt, unit, testflag)
      type(prt_t), intent(in) :: prt
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine prt_write
    elemental module function prt_match (prt1, prt2) result (match)
      logical :: match
      type(prt_t), intent(in) :: prt1, prt2
    end function prt_match
    module function are_disjoint (prt_in1, prt_in2) result (flag)
      logical :: flag
      type(prt_t), intent(in) :: prt_in1, prt_in2
    end function are_disjoint
    module subroutine subevt_init (subevt, n_active)
      type(subevt_t), intent(out) :: subevt
      integer, intent(in), optional :: n_active
    end subroutine subevt_init
    module subroutine subevt_reset (subevt, n_active)
      class(subevt_t), intent(inout) :: subevt
      integer, intent(in) :: n_active
    end subroutine subevt_reset
    module subroutine subevt_write (object, unit, prefix, pacified)
      class(subevt_t), intent(in) :: object
      integer, intent(in), optional :: unit
      character(*), intent(in), optional :: prefix
      logical, intent(in), optional :: pacified
    end subroutine subevt_write
    module subroutine subevt_assign (subevt, subevt_in)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: subevt_in
    end subroutine subevt_assign
    module subroutine subevt_set_beam (subevt, i, pdg, p, p2, src)
      class(subevt_t), intent(inout) :: subevt
      integer, intent(in) :: i
      integer, intent(in) :: pdg
      type(vector4_t), intent(in) :: p
      real(default), intent(in) :: p2
      integer, dimension(:), intent(in), optional :: src
    end subroutine subevt_set_beam
    module subroutine subevt_set_incoming (subevt, i, pdg, p, p2, src)
      class(subevt_t), intent(inout) :: subevt
      integer, intent(in) :: i
      integer, intent(in) :: pdg
      type(vector4_t), intent(in) :: p
      real(default), intent(in) :: p2
      integer, dimension(:), intent(in), optional :: src
    end subroutine subevt_set_incoming
    module subroutine subevt_set_outgoing (subevt, i, pdg, p, p2, src)
      class(subevt_t), intent(inout) :: subevt
      integer, intent(in) :: i
      integer, intent(in) :: pdg
      type(vector4_t), intent(in) :: p
      real(default), intent(in) :: p2
      integer, dimension(:), intent(in), optional :: src
    end subroutine subevt_set_outgoing
    module subroutine subevt_set_composite (subevt, i, p, src)
      class(subevt_t), intent(inout) :: subevt
      integer, intent(in) :: i
      type(vector4_t), intent(in) :: p
      integer, dimension(:), intent(in) :: src
    end subroutine subevt_set_composite
    module subroutine subevt_set_pdg_beam (subevt, pdg)
      class(subevt_t), intent(inout) :: subevt
      integer, dimension(:), intent(in) :: pdg
    end subroutine subevt_set_pdg_beam
    module subroutine subevt_set_pdg_incoming (subevt, pdg)
      class(subevt_t), intent(inout) :: subevt
      integer, dimension(:), intent(in) :: pdg
    end subroutine subevt_set_pdg_incoming
    module subroutine subevt_set_pdg_outgoing (subevt, pdg)
      class(subevt_t), intent(inout) :: subevt
      integer, dimension(:), intent(in) :: pdg
    end subroutine subevt_set_pdg_outgoing
    module subroutine subevt_set_p_beam (subevt, p)
      class(subevt_t), intent(inout) :: subevt
      type(vector4_t), dimension(:), intent(in) :: p
    end subroutine subevt_set_p_beam
    module subroutine subevt_set_p_incoming (subevt, p)
      class(subevt_t), intent(inout) :: subevt
      type(vector4_t), dimension(:), intent(in) :: p
    end subroutine subevt_set_p_incoming
    module subroutine subevt_set_p_outgoing (subevt, p)
      class(subevt_t), intent(inout) :: subevt
      type(vector4_t), dimension(:), intent(in) :: p
    end subroutine subevt_set_p_outgoing
    module subroutine subevt_set_p2_beam (subevt, p2)
      class(subevt_t), intent(inout) :: subevt
      real(default), dimension(:), intent(in) :: p2
    end subroutine subevt_set_p2_beam
    module subroutine subevt_set_p2_incoming (subevt, p2)
      class(subevt_t), intent(inout) :: subevt
      real(default), dimension(:), intent(in) :: p2
    end subroutine subevt_set_p2_incoming
    module subroutine subevt_set_p2_outgoing (subevt, p2)
      class(subevt_t), intent(inout) :: subevt
      real(default), dimension(:), intent(in) :: p2
    end subroutine subevt_set_p2_outgoing
    module subroutine subevt_polarize (subevt, i, h)
      type(subevt_t), intent(inout) :: subevt
      integer, intent(in) :: i, h
    end subroutine subevt_polarize
    module subroutine subevt_colorize (subevt, i, col, acl)
      type(subevt_t), intent(inout) :: subevt
      integer, intent(in) :: i, col, acl
    end subroutine subevt_colorize
    module function subevt_is_nonempty (subevt) result (flag)
      logical :: flag
      class(subevt_t), intent(in) :: subevt
    end function subevt_is_nonempty
    module function subevt_get_length (subevt) result (length)
      integer :: length
      class(subevt_t), intent(in) :: subevt
    end function subevt_get_length
    module function subevt_get_prt (subevt, i) result (prt)
      type(prt_t) :: prt
      class(subevt_t), intent(in) :: subevt
      integer, intent(in) :: i
    end function subevt_get_prt
    module function subevt_get_sqrts_hat (subevt) result (sqrts_hat)
      class(subevt_t), intent(in) :: subevt
      real(default) :: sqrts_hat
    end function subevt_get_sqrts_hat
    module function subevt_get_n_in (subevt) result (n_in)
      class(subevt_t), intent(in) :: subevt
      integer :: n_in
    end function subevt_get_n_in
    module function subevt_get_n_out (subevt) result (n_out)
      class(subevt_t), intent(in) :: subevt
      integer :: n_out
    end function subevt_get_n_out
    module function c_prt_from_subevt (subevt, i) result (c_prt)
      type(c_prt_t) :: c_prt
      type(subevt_t), intent(in) :: subevt
      integer, intent(in) :: i
    end function c_prt_from_subevt
    module function c_prt_array_from_subevt (subevt) result (c_prt_array)
      type(subevt_t), intent(in) :: subevt
      type(c_prt_t), dimension(subevt%n_active) :: c_prt_array
    end function c_prt_array_from_subevt
    module subroutine subevt_join (subevt, pl1, pl2, mask2)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl1, pl2
      logical, dimension(:), intent(in), optional :: mask2
    end subroutine subevt_join
    module subroutine subevt_combine (subevt, pl1, pl2, mask12)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl1, pl2
      logical, dimension(:,:), intent(in), optional :: mask12
    end subroutine subevt_combine
    module subroutine subevt_collect (subevt, pl1, mask1)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl1
      logical, dimension(:), intent(in) :: mask1
    end subroutine subevt_collect
    module subroutine subevt_cluster (subevt, pl1, dcut, mask1, jet_def, &
         keep_jets, exclusive)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl1
      real(default), intent(in) :: dcut
      logical, dimension(:), intent(in) :: mask1
      type(jet_definition_t), intent(in) :: jet_def
      logical, intent(in) :: keep_jets, exclusive
    end subroutine subevt_cluster
    module subroutine subevt_recombine (subevt, pl, mask1, reco_r0, keep_flv)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl
      logical, dimension(:), intent(in) :: mask1
      logical, intent(in) :: keep_flv
      real(default), intent(in) :: reco_r0
    end subroutine subevt_recombine
    module subroutine subevt_select (subevt, pl, mask1)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl
      logical, dimension(:), intent(in) :: mask1
    end subroutine subevt_select
    module subroutine subevt_extract (subevt, pl, index)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl
      integer, intent(in) :: index
    end subroutine subevt_extract
    module subroutine subevt_sort_pdg (subevt, pl)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl
    end subroutine subevt_sort_pdg
    module subroutine subevt_sort_int (subevt, pl, ival)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl
      integer, dimension(:), intent(in) :: ival
    end subroutine subevt_sort_int
    module subroutine subevt_sort_real (subevt, pl, rval)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl
      real(default), dimension(:), intent(in) :: rval
    end subroutine subevt_sort_real
    module subroutine subevt_select_pdg_code (subevt, aval, subevt_in, prt_type)
      type(subevt_t), intent(inout) :: subevt
      type(pdg_array_t), intent(in) :: aval
      type(subevt_t), intent(in) :: subevt_in
      integer, intent(in), optional :: prt_type
    end subroutine subevt_select_pdg_code
    module subroutine pacify_prt (prt)
      class(prt_t), intent(inout) :: prt
    end subroutine pacify_prt
    module subroutine pacify_subevt (subevt)
      class(subevt_t), intent(inout) :: subevt
    end subroutine pacify_subevt
  end interface

end module subevents
