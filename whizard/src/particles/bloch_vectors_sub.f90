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

submodule (bloch_vectors) bloch_vectors_s

  use physics_defs, only: SCALAR, SPINOR, VECTOR, VECTORSPINOR, TENSOR
  use su_algebra

  implicit none

contains

  function bloch_factor (s) result (f)
    real(default) :: f
    integer, intent(in) :: s
    select case (s)
    case (SCALAR)
       f = 0
    case (SPINOR)
       f = 1
    case (VECTOR)
       f = 2 * sqrt (3._default) / 3
    case (VECTORSPINOR)
       f = 2 * sqrt (6._default) / 4
    case (TENSOR)
       f = 2 * sqrt (10._default) / 5
    case default
       f = 0
    end select
  end function bloch_factor

  module subroutine bloch_vector_init_unpolarized (pol, spin_type)
    class(bloch_vector_t), intent(out) :: pol
    integer, intent(in) :: spin_type
    pol%spin_type = spin_type
  end subroutine bloch_vector_init_unpolarized

  module subroutine bloch_vector_init (pol, spin_type)
    class(bloch_vector_t), intent(out) :: pol
    integer, intent(in) :: spin_type
    pol%spin_type = spin_type
    select case (spin_type)
    case (SCALAR,SPINOR,VECTOR,VECTORSPINOR,TENSOR)
       allocate (pol%a (algebra_dimension (spin_type)), source = 0._default)
    end select
  end subroutine bloch_vector_init

  module subroutine bloch_vector_from_array (pol, a)
    class(bloch_vector_t), intent(inout) :: pol
    real(default), dimension(:), allocatable, intent(in) :: a
    pol%a(:) = a
  end subroutine bloch_vector_from_array

  module subroutine bloch_vector_to_array (pol, a)
    class(bloch_vector_t), intent(in) :: pol
    real(default), dimension(:), allocatable, intent(out) :: a
    if (pol%is_defined ())  allocate (a (size (pol%a)), source = pol%a)
  end subroutine bloch_vector_to_array

  module subroutine bloch_vector_write_raw (pol, u)
    class(bloch_vector_t), intent(in) :: pol
    integer, intent(in) :: u
    write (u) pol%spin_type
    write (u) allocated (pol%a)
    if (allocated (pol%a)) then
       write (u) pol%a
    end if
  end subroutine bloch_vector_write_raw

  module subroutine bloch_vector_read_raw (pol, u, iostat)
    class(bloch_vector_t), intent(out) :: pol
    integer, intent(in) :: u
    integer, intent(out) :: iostat
    integer :: s
    logical :: polarized
    read (u, iostat=iostat) s
    read (u, iostat=iostat) polarized
    if (iostat /= 0)  return
    if (polarized) then
       call pol%init (s)
       read (u, iostat=iostat) pol%a
    else
       call pol%init_unpolarized (s)
    end if
  end subroutine bloch_vector_read_raw

  module function get_n_states (pol) result (n)
    class(bloch_vector_t), intent(in) :: pol
    integer :: n
    n = fundamental_dimension (pol%spin_type)
  end function get_n_states

  module function get_length (pol) result (n)
    class(bloch_vector_t), intent(in) :: pol
    integer :: n
    n = algebra_dimension (pol%spin_type)
  end function get_length

  module function bv_helicity_index (pol, h) result (i)
    class(bloch_vector_t), intent(in) :: pol
    integer, intent(in) :: h
    integer :: i
    i = helicity_index (pol%spin_type, h)
  end function bv_helicity_index

  module function bv_helicity_value (pol, i) result (h)
    class(bloch_vector_t), intent(in) :: pol
    integer, intent(in) :: i
    integer :: h
    h = helicity_value (pol%spin_type, i)
  end function bv_helicity_value

  module function bv_factor (pol) result (f)
    class(bloch_vector_t), intent(in) :: pol
    real(default) :: f
    f = bloch_factor (pol%spin_type)
  end function bv_factor

  module function bloch_vector_is_defined (pol) result (flag)
    class(bloch_vector_t), intent(in) :: pol
    logical :: flag
    flag = pol%spin_type /= UNKNOWN
  end function bloch_vector_is_defined

  module function bloch_vector_is_polarized (pol) result (flag)
    class(bloch_vector_t), intent(in) :: pol
    logical :: flag
    flag = allocated (pol%a)
  end function bloch_vector_is_polarized

  module function bloch_vector_is_diagonal (pol) result (diagonal)
    class(bloch_vector_t), intent(in) :: pol
    logical :: diagonal
    integer :: s, i
    s = pol%spin_type
    diagonal = .true.
    if (pol%is_polarized ()) then
       do i = 1, size (pol%a)
          if (is_cartan_generator (s, i))  cycle
          if (pol%a(i) /= 0) then
             diagonal = .false.
             return
          end if
       end do
    end if
  end function bloch_vector_is_diagonal

  module function bloch_vector_get_norm (pol) result (norm)
    class(bloch_vector_t), intent(in) :: pol
    real(default) :: norm
    select case (pol%spin_type)
    case (SPINOR,VECTOR,VECTORSPINOR,TENSOR)
       norm = sqrt (dot_product (pol%a, pol%a))
    case default
       norm = 1
    end select
  end function bloch_vector_get_norm

  module subroutine bloch_vector_init_diagonal (pol, spin_type, rd)
    class(bloch_vector_t), intent(out) :: pol
    integer, intent(in) :: spin_type
    real(default), dimension(:), intent(in) :: rd
    call pol%init (spin_type)
    call pol%set (rd)
  end subroutine bloch_vector_init_diagonal

  module subroutine bloch_vector_set_diagonal (pol, rd)
    class(bloch_vector_t), intent(inout) :: pol
    real(default), dimension(:), intent(in) :: rd
    integer :: s
    s = pol%spin_type
    select case (s)
    case (SCALAR,SPINOR,VECTOR,VECTORSPINOR,TENSOR)
       pol%a(:) = cartan_coeff (s, rd) / bloch_factor (s)
    end select
  end subroutine bloch_vector_set_diagonal

  module subroutine bloch_vector_init_max_weight (pol, spin_type)
    class(bloch_vector_t), intent(out) :: pol
    integer, intent(in) :: spin_type
    call pol%init (spin_type)
    select case (spin_type)
    case (VECTOR)
       call pol%set ([0.5_default, 0._default, 0.5_default])
    case (VECTORSPINOR)
       call pol%set ([0.5_default, 0._default, 0._default, 0.5_default])
    case (TENSOR)
       call pol%set ([0.5_default, 0._default, 0._default, 0._default, 0.5_default])
    end select
  end subroutine bloch_vector_init_max_weight

  module subroutine bloch_vector_init_vector (pol, s, a)
    class(bloch_vector_t), intent(out) :: pol
    integer, intent(in) :: s
    real(default), dimension(3), intent(in) :: a
    call pol%init_max_weight (s)
    select case (s)
    case (SPINOR, VECTOR, VECTORSPINOR, TENSOR)
       pol%a(1:3) = a / bloch_factor (s)
    end select
  end subroutine bloch_vector_init_vector

  module subroutine bloch_vector_to_vector (pol, a)
    class(bloch_vector_t), intent(in) :: pol
    real(default), dimension(3), intent(out) :: a
    integer :: s
    s = pol%spin_type
    select case (s)
    case (SPINOR, VECTOR, VECTORSPINOR, TENSOR)
       a = pol%a(1:3) * bloch_factor (s)
    case default
       a = 0
    end select
  end subroutine bloch_vector_to_vector

  module subroutine bloch_vector_init_matrix (pol, spin_type, r)
    class(bloch_vector_t), intent(out) :: pol
    integer, intent(in) :: spin_type
    complex(default), dimension(:,:), intent(in) :: r
    select case (spin_type)
    case (SCALAR,SPINOR,VECTOR,VECTORSPINOR,TENSOR)
       call pol%init (spin_type)
       call pol%set (r)
    case default
       call pol%init (UNKNOWN)
    end select
  end subroutine bloch_vector_init_matrix

  module subroutine bloch_vector_set_matrix (pol, r)
    class(bloch_vector_t), intent(inout) :: pol
    complex(default), dimension(:,:), intent(in) :: r
    real(default), dimension(:), allocatable :: rd
    integer :: s, d, i, j, h1, h2, ir, ii
    s = pol%spin_type
    select case (s)
    case (SCALAR,SPINOR,VECTOR,VECTORSPINOR,TENSOR)
       d = fundamental_dimension (s)
       allocate (rd (d))
       do i = 1, d
          rd(i) = r(i,i)
       end do
       call pol%set (rd)
       do i = 1, d
          h1 = helicity_value (s, i)
          do j = i+1, d
             h2 = helicity_value (s, j)
             ir = root_index (s, h1, h2, .true.)
             ii = root_index (s, h1, h2, .false.)
             pol%a(ir) = real (r(j,i) + r(i,j)) / bloch_factor (s)
             pol%a(ii) = aimag (r(j,i) - r(i,j)) / bloch_factor (s)
          end do
       end do
    end select
  end subroutine bloch_vector_set_matrix

  module subroutine bloch_vector_to_matrix (pol, r, only_max_weight)
    class(bloch_vector_t), intent(in) :: pol
    complex(default), dimension(:,:), intent(out), allocatable :: r
    logical, intent(in), optional :: only_max_weight
    integer :: d, s, h0, ng, ai, h, h1, h2, i, j
    logical :: is_real, only_max
    complex(default) :: val
    if (.not. pol%is_polarized ())  return
    s = pol%spin_type
    only_max = .false.
    select case (s)
    case (VECTOR, VECTORSPINOR, TENSOR)
       if (present (only_max_weight))  only_max = only_max_weight
    end select
    if (only_max) then
       ng = 2
       h0 = helicity_value (s, 1)
    else
       ng = algebra_dimension (s)
       h0 = 0
    end if
    d = fundamental_dimension (s)
    allocate (r (d, d), source = (0._default, 0._default))
    do i = 1, d
       h = helicity_value (s, i)
       if (abs (h) < h0)  cycle
       r(i,i) = 1._default / d &
            + dot_product (cartan_element (s, h), pol%a) * bloch_factor (s)
    end do
    do ai = 1, ng
       if (is_cartan_generator (s, ai))  cycle
       call root_helicity (s, ai, h1, h2, is_real)
       i = helicity_index (s, h1)
       j = helicity_index (s, h2)
       if (is_real) then
          val = cmplx (pol%a(ai) / 2 * bloch_factor (s), 0._default, &
               kind=default)
          r(i,j) = r(i,j) + val
          r(j,i) = r(j,i) + val
       else
          val = cmplx (0._default, pol%a(ai) / 2 * bloch_factor (s), &
               kind=default)
          r(i,j) = r(i,j) - val
          r(j,i) = r(j,i) + val
       end if
    end do
  end subroutine bloch_vector_to_matrix


end submodule bloch_vectors_s

