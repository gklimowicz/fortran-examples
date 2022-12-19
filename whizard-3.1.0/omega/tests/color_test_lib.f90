! color_test_lib.f90 --
! color_test_lib.f90 -- O'Mega self test support
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     Christian Speckner <cnspeckn@googlemail.com>
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

module color_test_lib

  use kinds
  use omega_color

  implicit none
  private

  real(kind=default), parameter :: TOLERANCE = 1

  public :: compare_color_flows

  interface match_color_flow
     ! module procedure match_color_flow_with_allocation
     module procedure match_color_flow_without_allocation
  end interface

contains

  subroutine compare_color_flows (process, &
       cflows_expected, gflags_expected, cfactors_expected, &
       p_prt_in, p_prt_out, n_cindices, n_cflows, color_flows, &
       n_cfactors, color_factors, ok)
    character(*), intent(in) :: process
    integer, dimension(:,:,:), intent(in) :: cflows_expected
    logical, dimension(:,:), intent(in) :: gflags_expected
    type(omega_color_factor), dimension(:), intent(in) :: cfactors_expected
    logical, intent(out) :: ok
    interface
       pure function p_prt_in () result (n)
         integer :: n
       end function p_prt_in
       pure function p_prt_out () result (n)
         integer :: n
       end function p_prt_out
       pure function n_cindices () result (n)
         integer :: n
       end function n_cindices
       pure function n_cflows () result (n)
         integer :: n
       end function n_cflows
       pure subroutine color_flows (a, g)
         integer, dimension(:,:,:), intent(out) :: a
         logical, dimension(:,:), intent(out) :: g
       end subroutine color_flows
       pure function n_cfactors () result (n)
         integer :: n
       end function n_cfactors
       pure subroutine color_factors (cf)
         use omega_color
         type(omega_color_factor), dimension(:), intent(out) :: cf
       end subroutine color_factors
    end interface

    integer, dimension(n_cindices (),p_prt_in()+p_prt_out(), n_cflows()) :: cflows
    logical, dimension(p_prt_in()+p_prt_out(), n_cflows()) :: gflags
    type(omega_color_factor), dimension(n_cfactors()) :: cfactors
    real(kind=default), dimension(n_cflows(), n_cflows()) :: &
         cfactor_table, cfactor_table_expected
    integer, dimension(size (cflows_expected,dim=3)) :: indices, identity
    integer :: i

    ok = .true.
    identity = (/ (i, i = 1, size (identity)) /)

    if (n_cindices () .ne. size (cflows_expected, dim = 1)) then
       print *, "MISMATCH #INDICES: ", n_cindices (), &
            " <> ", size (cflows_expected, dim = 1)
       ok = .false.
    end if
    if (p_prt_in () + p_prt_out () .ne. size (cflows_expected, dim = 2)) then
       print *, "MISMATCH #PARTICLES: ", p_prt_in () + p_prt_out (), &
            " <> ", size (cflows_expected, dim = 2)
       ok = .false.
    end if
    if (n_cflows () .ne. size (cflows_expected, dim = 3)) then
       print *, "MISMATCH #COLORFLOWS: ", n_cflows (), &
            " <> ", size (cflows_expected, dim = 3)
       ok = .false.
    end if
    if (n_cfactors () .ne. size (cfactors_expected)) then
       print *, "MISMATCH #COLORFACTORS: ", n_cfactors (), &
            " <> ", size (cfactors_expected)
       ok = .false.
    end if

    call color_flows (cflows, gflags)
    call match_color_flows (cflows_expected, gflags_expected, cflows, gflags, indices)
    if (any (indices .le. 0)) then
       print *, "COLOR FLOW MISMATCH: "
       print *, "indices = ", indices
       ok = .false.
    end if

    call color_factors (cfactors)
    call unpack_color_factors (cfactors, indices, cfactor_table)
    call unpack_color_factors (cfactors_expected, identity, cfactor_table_expected)
    if (      any (abs (cfactor_table - cfactor_table_expected) &
         .gt. TOLERANCE * epsilon (real (1, kind=default)))) then
       print *, "COLOR FACTOR MISMATCH: "
       print *, "expected = ", cfactor_table_expected
       print *, "computed = ", cfactor_table
       ok = .false.
    end if

    if (ok) then
       print *, "color flow/factor test for process ", process, " passed."
    else
       print *, "color flow/factor test for process ", process, " failed!"
    end if

  end subroutine compare_color_flows

  subroutine unpack_color_factors (cfactors, indices, a)
    type(omega_color_factor), dimension(:), intent(in) :: cfactors
    integer, dimension(:), intent(in) :: indices
    real(kind=default), dimension(:,:), intent(out) :: a
    integer :: i, j
    if (      (size (a, dim = 1) .ne. size (a, dim = 2)) &
         .or. (size (a, dim = 1) .ne. size (indices))) then
       print *, "unpack_color_factors: #mismatch"
       stop 1
    end if
    a = 0
    ! This could be made more efficient by reverting INDICES first
    do i = 1, size (indices)
       do j = 1, size (indices)
          a(i,j) = find_color_factor (cfactors, indices(i), indices(j))
       end do
    end do
  end subroutine unpack_color_factors

  function find_color_factor (cfactors, i1, i2) result (factor)
    real(kind=default) :: factor
    type(omega_color_factor), dimension(:), intent(in) :: cfactors
    integer, intent(in) :: i1, i2
    integer :: i
    factor = 0
    do i = 1, size (cfactors)
       if ((cfactors(i)%i1 .eq. i1) .and. (cfactors(i)%i2 .eq. i2)) then
          factor = cfactors(i)%factor
          return
       end if
    end do
  end function find_color_factor

  subroutine match_color_flows (cfs1, gfs1, cfs2, gfs2, indices)
    integer, dimension(:,:,:), intent(in) :: cfs1, cfs2
    logical, dimension(:,:), intent(in) :: gfs1, gfs2
    integer, dimension(:), intent(out) :: indices
    integer :: i
    if (      (size (cfs1, dim = 3) .ne. size (gfs1, dim = 2)) &
         .or. (size (cfs2, dim = 3) .ne. size (gfs2, dim = 2)) &
         .or. (size (cfs1, dim = 3) .ne. size (indices))) then
       print *, "match_color_flows: length mismatch"
       stop 1
    end if
    indices = -1
    if (size (cfs1, dim = 3) .eq. size (cfs2, dim = 3)) then
       do i = 1, size (cfs1, dim = 3)
          indices(i) = find_color_flow (cfs1(:,:,i), gfs1(:,i), cfs2, gfs2)
       end do
    end if
  end subroutine match_color_flows

  function find_color_flow (cf, gf, cfs, gfs) result (n)
    integer :: n
    integer, dimension(:,:), intent(in) :: cf
    logical, dimension(:), intent(in) :: gf
    integer, dimension(:,:,:), intent(in) :: cfs
    logical, dimension(:,:), intent(in) :: gfs
    integer :: i
    if (size (cfs, dim = 3) .ne. size (gfs, dim = 2)) then
       print *, "find_color_flows: #mismatch"
       stop 1
    end if
    do i = 1, size (cfs, dim = 3)
       if (match_color_flow (cf, gf, cfs(:,:,i), gfs(:,i))) then
          n = i
          return
       end if
    end do
    n = -1
  end function find_color_flow

  function match_color_flow_without_allocation (cf1, gf1, cf2, gf2) result (match)
    logical :: match
    integer, dimension(:,:), intent(in) :: cf1, cf2
    logical, dimension(:), intent(in) :: gf1, gf2
    integer :: n, i, j
    match = .false.
    if (      (size (cf1, dim = 1) .ne. 2) &
         .or. (size (cf2, dim = 1) .ne. 2)) then
       print *, "match_color_flow: wrong #indices"
       stop 1
    end if
    if (      (size (cf1, dim = 2) .ne. size (gf1)) &
         .or. (size (cf2, dim = 2) .ne. size (gf2)) &
         .or. (size (cf1, dim = 2) .ne. size (cf2, dim = 2))) then
       print *, "match_color_flow: length mismatch"
       stop 1
    end if
    if (all (gf1 .eqv. gf2)) then
       n = size (cf1, dim = 2)
       match = .true.
       do i = 1, n
          do j = 1, n
             ! compare incidence matrix elements and immediately
             ! return .false., if one doesn't match
             if (        incidence_matrix_element (cf1, i, j) &
                  .neqv. incidence_matrix_element (cf2, i, j)) then
                match = .false.
                return
             end if
          end do
       end do
    else
       match = .false.
    end if
  end function match_color_flow_without_allocation

  pure function incidence_matrix_element (cf, i, j) result (linked)
    integer, dimension(:,:), intent(in) :: cf
    integer, intent(in) :: i, j
    logical :: linked
    linked = (i .ne. j) .and. &
      (      (cf(1,i) .ne. 0) .and. (cf(1,i) .eq. -cf(2,j)) &
        .or. (cf(1,i) .ne. 0) .and. (cf(1,i) .eq.  cf(1,j)) &
        .or. (cf(2,i) .ne. 0) .and. (cf(2,i) .eq.  cf(2,j)))
  end function incidence_matrix_element

  function match_color_flow_with_allocation (cf1, gf1, cf2, gf2) result (match)
    logical :: match
    integer, dimension(:,:), intent(in) :: cf1, cf2
    logical, dimension(:), intent(in) :: gf1, gf2
    logical, dimension(size(cf1,dim=2),size(cf1,dim=2)) :: a1, a2
    integer :: n
    match = .false.
    if (      (size (cf1, dim = 1) .ne. 2) &
         .or. (size (cf2, dim = 1) .ne. 2)) then
       print *, "match_color_flow: wrong #indices"
       stop 1
    end if
    if (      (size (cf1, dim = 2) .ne. size (gf1)) &
         .or. (size (cf2, dim = 2) .ne. size (gf2)) &
         .or. (size (cf1, dim = 2) .ne. size (cf2, dim = 2))) then
       print *, "match_color_flow: length mismatch"
       stop 1
    end if
    if (all (gf1 .eqv. gf2)) then
       n = size(cf1,dim=2)
       call incidence_matrix_of_color_flow (a1, cf1)
       call incidence_matrix_of_color_flow (a2, cf2)
       ! print *, "a1: ", a1
       ! print *, "a2: ", a2
       if (all (a1 .eqv. a2)) then
          match = .true.
       end if
    end if
  end function match_color_flow_with_allocation

  pure subroutine incidence_matrix_of_color_flow (a, cf)
    logical, dimension(:,:), intent (out) :: a
    integer, dimension(:,:), intent(in) :: cf
    integer :: n, i, j
    n = size (cf, dim = 2)
    forall (i = 1:n, j = 1:n)
       a(i,j) = incidence_matrix_element (cf, i, j)
    end forall
  end subroutine incidence_matrix_of_color_flow

end module color_test_lib
