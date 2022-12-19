!$Id: fastjet.f90 6133 2014-09-17 14:42:33Z kilian $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

submodule (fastjet) fastjet_s

  implicit none

contains

  ! Fastjet banner, print explicitly to control order of output
  module subroutine print_banner ()
    call fastjet_print_banner ()
  end subroutine print_banner

  ! Procedures for pseudojets
  module subroutine pseudojet_init (j, px, py, pz, E)
    class(pseudojet_t), intent(out) :: j
    real(default), intent(in) :: px, py, pz, E
    real(c_double) :: jx = 0
    real(c_double) :: jy = 0
    real(c_double) :: jz = 0
    real(c_double) :: je = 0
    jx = px
    jy = py
    jz = pz
    je = E
    j%cptr = new_pseudojet (jx, jy, jz, je)
  end subroutine pseudojet_init

  module subroutine pseudojet_init_vector4 (j, p)
    class(pseudojet_t), intent(out) :: j
    type(vector4_t), intent(in) :: p
    real(c_double) :: jx = 0
    real(c_double) :: jy = 0
    real(c_double) :: jz = 0
    real(c_double) :: je = 0
    jx = vector4_get_component (p, 1)
    jy = vector4_get_component (p, 2)
    jz = vector4_get_component (p, 3)
    je = vector4_get_component (p, 0)
    j%cptr = new_pseudojet (jx, jy, jz, je)
  end subroutine pseudojet_init_vector4

  module subroutine pseudojet_final (j)
    class(pseudojet_t), intent(inout) :: j
    call pseudojet_delete (j%cptr)
    j%cptr = c_null_ptr
  end subroutine pseudojet_final

  module function pseudojet_e (j) result (e)
    class(pseudojet_t), intent(in) :: j
    real(default) :: e
    e = pseudojet_get_e (j%cptr)
  end function pseudojet_e

  module function pseudojet_px (j) result (p)
    class(pseudojet_t), intent(in) :: j
    real(default) :: p
    p = pseudojet_get_px (j%cptr)
  end function pseudojet_px

  module function pseudojet_py (j) result (p)
    class(pseudojet_t), intent(in) :: j
    real(default) :: p
    p = pseudojet_get_py (j%cptr)
  end function pseudojet_py

  module function pseudojet_pz (j) result (p)
    class(pseudojet_t), intent(in) :: j
    real(default) :: p
    p = pseudojet_get_pz (j%cptr)
  end function pseudojet_pz

  module function pseudojet_perp (j) result (p)
    class(pseudojet_t), intent(in) :: j
    real(default) :: p
    p = pseudojet_get_perp (j%cptr)
  end function pseudojet_perp

  module function pseudojet_rap (j) result (p)
    class(pseudojet_t), intent(in) :: j
    real(default) :: p
    p = pseudojet_get_rap (j%cptr)
  end function pseudojet_rap

  module function pseudojet_phi (j) result (p)
    class(pseudojet_t), intent(in) :: j
    real(default) :: p
    p = pseudojet_get_phi (j%cptr)
  end function pseudojet_phi

  module function pseudojet_constituents (j) result (prt)
    class(pseudojet_t), intent(in) :: j
    type(pseudojet_vector_t) :: prt
    prt%cptr = pseudojet_get_constituents (j%cptr)
  end function

  module function pseudojet_contains_prt (j, prt) result (flag)
    class(pseudojet_t), intent(in) :: j
    type(pseudojet_t), intent(in) :: prt
    logical :: flag
    flag = pseudojet_contains (j%cptr, prt%cptr)
  end function pseudojet_contains_prt


  ! Procedures for pseudojet vectors
  module function pseudojet_vector (j) result (jv)
    type(pseudojet_t), dimension(:), intent(in) :: j
    type(pseudojet_vector_t) :: jv
    call jv%init (j)
  end function pseudojet_vector

  module subroutine pseudojet_vector_init (jv, j)
    class(pseudojet_vector_t), intent(out) :: jv
    type(pseudojet_t), dimension(:), intent(in) :: j
    type(c_ptr), dimension(:), allocatable :: cptr
    allocate (cptr (size (j)))
    cptr = j%cptr
    jv%cptr = new_pseudojet_vector (cptr, size (j))
  end subroutine pseudojet_vector_init

  module subroutine pseudojet_vector_final (jv)
    class(pseudojet_vector_t), intent(inout) :: jv
    call pseudojet_vector_delete (jv%cptr)
    jv%cptr = c_null_ptr
  end subroutine pseudojet_vector_final

  module function pseudojet_vector_size (jv) result (n)
    class(pseudojet_vector_t), intent(in) :: jv
    integer :: n
    n = pseudojet_vector_get_size (jv%cptr)
  end function pseudojet_vector_size

  module function pseudojet_vector_get (jv, i) result (j)
    class(pseudojet_vector_t), intent(in) :: jv
    integer, intent(in) :: i
    type(pseudojet_t) :: j
    j%cptr = pseudojet_vector_get_jet (jv%cptr, int(i-1, c_int))
  end function pseudojet_vector_get

  module subroutine pseudojet_array_from_vector (j, jv)
    type(pseudojet_t), dimension(:), allocatable, intent(out) :: j
    type(pseudojet_vector_t), intent(in) :: jv
    integer :: i, n
    n = jv%size ()
    allocate (j (n))
    do i = 1, n
       j(i) = jv%get (i)
    end do
  end subroutine pseudojet_array_from_vector

  module function sorted_by_pt (jets) result (sorted_jets)
    class(pseudojet_vector_t), intent(in) :: jets
    type(pseudojet_vector_t) :: sorted_jets
    sorted_jets%cptr = pseudojet_vector_sorted_by_pt (jets%cptr)
  end function sorted_by_pt

  ! Procedures for jet definitions
  module subroutine jet_definition_init (jet_def, jet_alg, r, p, jet_ycut)
    class(jet_definition_t), intent(out) :: jet_def
    integer(jet_algorithm_kind), intent(in) :: jet_alg
    real(default), intent(in) :: r
    real(default), intent(in), optional :: p
    real(default), intent(in), optional :: jet_ycut
    type(cpp_string_t) :: description_str
    real(default) :: ycut
    real(default) :: pp
    ycut = -1._default
    pp = -1._default
    if (present (jet_ycut)) then
       ycut = jet_ycut
    end if
    if (present (p)) then
       pp = p
    end if
    jet_def%cptr = new_jet_definition (jet_alg, real (r, c_double), &
         real (pp, c_double), real (ycut, c_double))
    call jet_def%description_str%init &
         (jet_definition_get_description (jet_def%cptr))
    jet_def%description_strlen = len (jet_def%description_str)
  end subroutine jet_definition_init

  module subroutine jet_definition_final (jet_def)
    class(jet_definition_t), intent(inout) :: jet_def
    call jet_def%description_str%final ()
    jet_def%description_strlen = 0
    call jet_definition_delete (jet_def%cptr)
    jet_def%cptr = c_null_ptr
  end subroutine jet_definition_final

  module function jet_definition_description (jet_def) result (string)
    class(jet_definition_t), intent(in) :: jet_def
    character(len=jet_def%description_strlen) :: string
    string = char (jet_def%description_str)
  end function jet_definition_description


  ! Procedures for cluster sequences
  module subroutine cluster_sequence_init (cs, jv, jet_def)
    class(cluster_sequence_t), intent(out) :: cs
    type(pseudojet_vector_t), intent(in) :: jv
    type(jet_definition_t), intent(in) :: jet_def
    cs%cptr = new_cluster_sequence (jv%cptr, jet_def%cptr)
  end subroutine cluster_sequence_init

  module subroutine cluster_sequence_final (cs)
    class(cluster_sequence_t), intent(inout) :: cs
    call cluster_sequence_delete (cs%cptr)
    cs%cptr = c_null_ptr
  end subroutine cluster_sequence_final

  module function cluster_sequence_inclusive_jets (cs) result (jets)
    class(cluster_sequence_t), intent(in) :: cs
    type(pseudojet_vector_t) :: jets
    jets%cptr = cluster_sequence_get_inclusive_jets (cs%cptr)
  end function cluster_sequence_inclusive_jets

  module function cluster_sequence_exclusive_jets (cs, dcut) result (jets)
    class(cluster_sequence_t), intent(in) :: cs
    type(pseudojet_vector_t) :: jets
    real(default), intent(in) :: dcut
    jets%cptr = cluster_sequence_get_exclusive_jets (cs%cptr, &
         real (dcut, c_double))
  end function cluster_sequence_exclusive_jets

  module subroutine cluster_sequence_assign_jet_indices (cs, jets, idx)
    class(cluster_sequence_t), intent(in) :: cs
    type(pseudojet_vector_t), intent(in) :: jets
    integer, dimension(:), intent(out) :: idx
    type(c_ptr) :: idx_cptr
    integer :: i
    idx_cptr = cluster_sequence_get_jet_indices (cs%cptr, jets%cptr)
    do i = 1, size (idx)
       idx(i) = int_vector_get (idx_cptr, int(i-1, c_int)) + 1
    end do
    call int_vector_delete (idx_cptr)
  end subroutine cluster_sequence_assign_jet_indices

end submodule fastjet_s

