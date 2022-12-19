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

module bloch_vectors_uti

  use kinds, only: default
  use physics_defs, only: UNKNOWN, SCALAR, SPINOR, VECTOR, VECTORSPINOR, TENSOR
  use su_algebra, only: algebra_dimension, fundamental_dimension, helicity_value

  use bloch_vectors

  implicit none
  private

  public :: bloch_vectors_1
  public :: bloch_vectors_2
  public :: bloch_vectors_3
  public :: bloch_vectors_4
  public :: bloch_vectors_5
  public :: bloch_vectors_6
  public :: bloch_vectors_7

contains

  subroutine bloch_vectors_1 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: bloch_vectors_1"
    write (u, "(A)")  "*   Purpose: test Bloch-vector &
         &polarization implementation"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization (unpolarized)"

    write (u, "(A)")
    write (u, "(A)")  "* unknown"
    call bloch_init (UNKNOWN)

    write (u, "(A)")
    write (u, "(A)")  "* s = 0"
    call bloch_init (SCALAR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call bloch_init (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call bloch_init (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call bloch_init (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call bloch_init (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: bloch_vectors_1"

  contains

    subroutine bloch_init (s)
      integer, intent(in) :: s
      type(bloch_vector_t) :: pol
      real(default), dimension(:), allocatable :: a
      integer :: i
      write (u, *)
      write (u, "(1X,L1,L1)", advance="no") &
           pol%is_defined (), pol%is_polarized ()
      call pol%init_unpolarized (s)
      write (u, "(1X,L1,L1)", advance="no") &
           pol%is_defined (), pol%is_polarized ()
      call pol%init (s)
      write (u, "(1X,L1,L1)", advance="no") &
           pol%is_defined (), pol%is_polarized ()
      write (u, *)
      call pol%to_array (a)
      if (allocated (a)) then
         write (u, "(*(F7.4))")  a
         a(:) = [(real (mod (i, 10), kind=default), i = 1, size (a))]
         call pol%from_array (a)
         call pol%to_array (a)
         write (u, "(*(F7.4))")  a
      else
         write (u, *)
         write (u, *)
      end if
    end subroutine bloch_init

  end subroutine bloch_vectors_1

  subroutine bloch_vectors_2 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: bloch_vectors_2"
    write (u, "(A)")  "*   Purpose: test Bloch-vector &
         &polarization implementation"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization (polarized, diagonal): &
         &display vector and norm"
    write (u, "(A)")  "*   transform back"

    write (u, "(A)")
    write (u, "(A)")  "* s = 0"
    call bloch_diagonal (SCALAR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call bloch_diagonal (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call bloch_diagonal (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call bloch_diagonal (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call bloch_diagonal (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: bloch_vectors_2"

  contains

    subroutine bloch_diagonal (s)
      integer, intent(in) :: s
      type(bloch_vector_t) :: pol
      real(default), dimension(:), allocatable :: a
      real(default), dimension(:), allocatable :: rd
      complex(default), dimension(:,:), allocatable :: r
      integer :: i, j, d
      real(default) :: rj
      real, parameter :: tolerance = 1.E-14_default
      d = fundamental_dimension (s)
      do i = 1, d
         allocate (rd (d), source = 0._default)
         rd(i) = 1
         call pol%init (s, rd)
         call pol%to_array (a)
         write (u, *)
         write (u, "(A,1X,I2)")  "h:", helicity_value (s, i)
         write (u, 1, advance="no")  a
         write (u, "(1X,L1)")  pol%is_diagonal ()
         write (u, 1)  pol%get_norm ()
         call pol%to_matrix (r)
         do j = 1, d
            rj = real (r(j,j))
            if (abs (rj) < tolerance)  rj = 0
            write (u, 1, advance="no")  rj
         end do
         write (u, "(1X,L1)")  matrix_is_diagonal (r)
         deallocate (a, rd, r)
      end do
1     format (99(1X,F7.4,:))
    end subroutine bloch_diagonal

    function matrix_is_diagonal (r) result (diagonal)
      complex(default), dimension(:,:), intent(in) :: r
      logical :: diagonal
      integer :: i, j
      diagonal = .true.
      do j = 1, size (r, 2)
         do i = 1, size (r, 1)
            if (i == j)  cycle
            if (r(i,j) /= 0) then
               diagonal = .false.
               return
            end if
         end do
      end do
    end function matrix_is_diagonal

  end subroutine bloch_vectors_2

  subroutine bloch_vectors_3 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: bloch_vectors_3"
    write (u, "(A)")  "*   Purpose: test Bloch-vector &
         &polarization implementation"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization (pure polarized, arbitrary):"
    write (u, "(A)")  "*   input matrix, transform, display norm, transform back"

    write (u, "(A)")
    write (u, "(A)")  "* s = 0"
    call bloch_arbitrary (SCALAR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call bloch_arbitrary (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call bloch_arbitrary (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call bloch_arbitrary (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call bloch_arbitrary (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: bloch_vectors_3"

  contains

    subroutine bloch_arbitrary (s)
      integer, intent(in) :: s
      type(bloch_vector_t) :: pol
      complex(default), dimension(:,:), allocatable :: r
      integer :: d
      d = fundamental_dimension (s)
      write (u, *)
      call init_matrix (d, r)
      where (abs (aimag (r)) < 1.e-14_default)  &
           r = cmplx (real(r, kind=default), 0._default, kind=default)
      call write_matrix (d, r)
      call pol%init (s, r)
      write (u, *)
      write (u, 2)  pol%get_norm (), pol%is_diagonal ()
      write (u, *)
      call pol%to_matrix (r)
      call write_matrix (d, r)
2     format (1X,F7.4,1X,L1)
    end subroutine bloch_arbitrary

    subroutine init_matrix (d, r)
      integer, intent(in) :: d
      complex(default), dimension(:,:), allocatable, intent(out) :: r
      complex(default), dimension(:), allocatable :: a
      real(default) :: norm
      integer :: i, j
      allocate (a (d))
      norm = 0
      do i = 1, d
         a(i) = cmplx (2*i-1, 2*i, kind=default)
         norm = norm + conjg (a(i)) * a(i)
      end do
      a = a / sqrt (norm)
      allocate (r (d,d))
      do i = 1, d
         do j = 1, d
            r(i,j) = conjg (a(i)) * a(j)
         end do
      end do
    end subroutine init_matrix

    subroutine write_matrix (d, r)
      integer, intent(in) :: d
      complex(default), dimension(:,:), intent(in) :: r
      integer :: i, j
      do i = 1, d
         do j = 1, d
            write (u, 1, advance="no")  r(i,j)
         end do
         write (u, *)
      end do
1     format (99(1X,'(',F7.4,',',F7.4,')',:))
    end subroutine write_matrix

  end subroutine bloch_vectors_3

  subroutine bloch_vectors_4 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: bloch_vectors_4"
    write (u, "(A)")  "*   Purpose: test Bloch-vector &
         &polarization implementation"
    write (u, "(A)")

    write (u, "(A)")  "* Raw I/O"

    write (u, "(A)")
    write (u, "(A)")  "* s = 0"
    call bloch_io (SCALAR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call bloch_io (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call bloch_io (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call bloch_io (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call bloch_io (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: bloch_vectors_4"

  contains

    subroutine bloch_io (s)
      integer, intent(in) :: s
      type(bloch_vector_t) :: pol
      real(default), dimension(:), allocatable :: a
      integer :: n, i, utmp, iostat
      n = algebra_dimension (s)
      allocate (a (n))
      a(:) = [(real (mod (i, 10), kind=default), i = 1, size (a))]
      write (u, *)
      write (u, "(*(F7.4))")  a
      call pol%init (s)
      call pol%from_array (a)
      open (newunit = utmp, status = "scratch", action = "readwrite", &
           form = "unformatted")
      call pol%write_raw (utmp)
      rewind (utmp)
      call pol%read_raw (utmp, iostat=iostat)
      close (utmp)
      call pol%to_array (a)
      write (u, "(*(F7.4))")  a
    end subroutine bloch_io

  end subroutine bloch_vectors_4

  subroutine bloch_vectors_5 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: bloch_vectors_5"
    write (u, "(A)")  "*   Purpose: test Bloch-vector &
         &polarization implementation"
    write (u, "(A)")

    write (u, "(A)")  "* Massless states: equipartition"

    write (u, "(A)")
    write (u, "(A)")  "* s = 0"
    call bloch_massless_unpol (SCALAR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call bloch_massless_unpol (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call bloch_massless_unpol (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call bloch_massless_unpol (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call bloch_massless_unpol (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: bloch_vectors_5"

  contains

    subroutine bloch_massless_unpol (s)
      integer, intent(in) :: s
      type(bloch_vector_t) :: pol
      complex(default), dimension(:,:), allocatable :: r
      real(default), dimension(:), allocatable :: a
      integer :: d
      d = fundamental_dimension (s)
      call pol%init_max_weight (s)
      call pol%to_matrix (r, only_max_weight = .false.)
      write (u, *)
      where (abs (r) < 1.e-14_default)  r = 0
      call write_matrix (d, r)
      call pol%to_matrix (r, only_max_weight = .true.)
      write (u, *)
      call write_matrix (d, r)
    end subroutine bloch_massless_unpol

    subroutine write_matrix (d, r)
      integer, intent(in) :: d
      complex(default), dimension(:,:), intent(in) :: r
      integer :: i, j
      do i = 1, d
         do j = 1, d
            write (u, 1, advance="no")  r(i,j)
         end do
         write (u, *)
      end do
1     format (99(1X,'(',F7.4,',',F7.4,')',:))
    end subroutine write_matrix

  end subroutine bloch_vectors_5

  subroutine bloch_vectors_6 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: bloch_vectors_6"
    write (u, "(A)")  "*   Purpose: test Bloch-vector &
         &polarization implementation"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization (pure polarized massless, arbitrary):"
    write (u, "(A)")  "*   input matrix, transform, display norm, transform back"

    write (u, "(A)")
    write (u, "(A)")  "* s = 0"
    call bloch_massless (SCALAR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call bloch_massless (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call bloch_massless (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call bloch_massless (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call bloch_massless (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: bloch_vectors_6"

  contains

    subroutine bloch_massless (s)
      integer, intent(in) :: s
      type(bloch_vector_t) :: pol
      complex(default), dimension(:,:), allocatable :: r
      integer :: d
      d = fundamental_dimension (s)
      write (u, *)
      call init_matrix (d, r)
      where (abs (aimag (r)) < 1.e-14_default)  &
           r = cmplx (real(r, kind=default), 0._default, kind=default)
      call write_matrix (d, r)
      call pol%init (s, r)
      write (u, *)
      write (u, 2)  pol%get_norm (), pol%is_diagonal ()
      write (u, *)
      call pol%to_matrix (r, only_max_weight = .true.)
      call write_matrix (d, r)
2     format (1X,F7.4,1X,L1)
    end subroutine bloch_massless

    subroutine init_matrix (d, r)
      integer, intent(in) :: d
      complex(default), dimension(:,:), allocatable, intent(out) :: r
      complex(default), dimension(:), allocatable :: a
      real(default) :: norm
      integer :: i, j
      allocate (a (d), source = (0._default, 0._default))
      norm = 0
      do i = 1, d, max (d-1, 1)
         a(i) = cmplx (2*i-1, 2*i, kind=default)
         norm = norm + conjg (a(i)) * a(i)
      end do
      a = a / sqrt (norm)
      allocate (r (d,d), source = (0._default, 0._default))
      do i = 1, d, max (d-1, 1)
         do j = 1, d, max (d-1, 1)
            r(i,j) = conjg (a(i)) * a(j)
         end do
      end do
    end subroutine init_matrix

    subroutine write_matrix (d, r)
      integer, intent(in) :: d
      complex(default), dimension(:,:), intent(in) :: r
      integer :: i, j
      do i = 1, d
         do j = 1, d
            write (u, 1, advance="no")  r(i,j)
         end do
         write (u, *)
      end do
1     format (99(1X,'(',F7.4,',',F7.4,')',:))
    end subroutine write_matrix

  end subroutine bloch_vectors_6

  subroutine bloch_vectors_7 (u)
    integer, intent(in) :: u

    write (u, "(A)")  "* Test output: bloch_vectors_7"
    write (u, "(A)")  "*   Purpose: test Bloch-vector &
         &polarization implementation"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization &
         &(pure polarized massless, arbitrary Bloch vector):"
    write (u, "(A)")  "*   input vector, transform, display norm, &
         &transform back"

    write (u, "(A)")
    write (u, "(A)")  "* s = 0"
    call bloch_massless_vector (SCALAR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1/2"
    call bloch_massless_vector (SPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 1"
    call bloch_massless_vector (VECTOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 3/2"
    call bloch_massless_vector (VECTORSPINOR)

    write (u, "(A)")
    write (u, "(A)")  "* s = 2"
    call bloch_massless_vector (TENSOR)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: bloch_vectors_7"

  contains

    subroutine bloch_massless_vector (s)
      integer, intent(in) :: s
      type(bloch_vector_t) :: pol
      real(default), dimension(3) :: a
      complex(default), dimension(:,:), allocatable :: r
      write (u, *)
      a =  [1._default, 2._default, 4._default]
      a = a / sqrt (sum (a ** 2))
      write (u, 2)  a
      call pol%init_vector (s, a)
      write (u, 2)  pol%get_norm ()
      call pol%to_vector (a)
      write (u, 2)  a
      call pol%to_matrix (r, only_max_weight = .false.)
      write (u, *)
      where (abs (r) < 1.e-14_default)  r = 0
      call write_matrix (r)
      call pol%to_matrix (r, only_max_weight = .true.)
      write (u, *)
      call write_matrix (r)
2     format (99(1X,F7.4,:))
    end subroutine bloch_massless_vector

    subroutine write_matrix (r)
      complex(default), dimension(:,:), intent(in) :: r
      integer :: i, j
      do i = 1, size (r, 1)
         do j = 1, size (r, 2)
            write (u, 1, advance="no")  r(i,j)
         end do
         write (u, *)
      end do
1     format (99(1X,'(',F7.4,',',F7.4,')',:))
    end subroutine write_matrix

  end subroutine bloch_vectors_7


end module bloch_vectors_uti
