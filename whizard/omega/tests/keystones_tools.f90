! keystones_tools.f90 -- tools for fusion/keystone/vertex tests
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 2019-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
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

module keystones_tools
  ! use ieee_arithmetic
  use kinds
  use constants
  ! use tao_random_numbers
  use omega95
  use omega95_bispinors, only: bispinor

  implicit none
  private

  public :: make_random

  interface make_random
     module procedure &
          make_random_real, &
          make_random_real_vector, &
          make_random_real_array, &
          make_random_complex, &
          make_random_complex_vector, &
          make_random_complex_array, &
          make_random_momentum, &
          make_random_momentum_vector, &
          make_random_vector, &
          make_random_vector_vector, &
          make_random_tensor, &
          make_random_tensor_vector, &
          make_random_tensor2odd, &
          make_random_tensor2odd_vector, &
          make_random_spinor, &
          make_random_spinor_vector, &
          make_random_conjspinor, &
          make_random_conjspinor_vector, &
          make_random_bispinor, &
          make_random_bispinor_vector
  end interface make_random

contains

  subroutine make_random_real (x, range)
    real(kind=default), intent(inout) :: x
    real(kind=default), intent(in), optional :: range
    call random_number (x)
    x = 2*x - 1
    if (present (range)) then
       x = range * x
    end if
  end subroutine make_random_real

  subroutine make_random_real_vector (x, range)
    real(kind=default), dimension(:), intent(inout) :: x
    real(kind=default), intent(in), optional :: range
    call random_number (x)
    x = 2*x - 1
    if (present (range)) then
       x = range * x
    end if
  end subroutine make_random_real_vector

  subroutine make_random_real_array (x, range)
    real(kind=default), dimension(:,:), intent(inout) :: x
    real(kind=default), intent(in), optional :: range
    call random_number (x)
    x = 2*x - 1
    if (present (range)) then
       x = range * x
    end if
  end subroutine make_random_real_array

  subroutine make_random_complex (z, range)
    complex(kind=default), intent(inout) :: z
    real(kind=default), intent(in), optional :: range
    real(kind=default) :: x, y
    call make_random_real (x, range)
    call make_random_real (y, range)
    z = cmplx (x, y, kind=default)
  end subroutine make_random_complex

  subroutine make_random_complex_vector (z, range)
    complex(kind=default), dimension(:), intent(inout) :: z
    real(kind=default), intent(in), optional :: range
    real(kind=default), dimension(size(z)) :: x, y
    call make_random_real_vector (x, range)
    call make_random_real_vector (y, range)
    z = cmplx (x, y, kind=default)
  end subroutine make_random_complex_vector
    
  subroutine make_random_complex_array (z, range)
    complex(kind=default), dimension(:,:), intent(inout) :: z
    real(kind=default), intent(in), optional :: range
    real(kind=default), dimension(size(z, dim=1),size(z, dim=2)) :: x, y
    call make_random_real_array (x, range)
    call make_random_real_array (y, range)
    z = cmplx (x, y, kind=default)
  end subroutine make_random_complex_array
    
  subroutine make_random_momentum (p, range)
    type(momentum), intent(inout) :: p
    real(kind=default), intent(in), optional :: range
    call make_random_real (p%t, range)
    call make_random_real_vector (p%x, range)
  end subroutine make_random_momentum

  subroutine make_random_momentum_vector (p, range)
    type(momentum), dimension(:), intent(inout) :: p
    real(kind=default), intent(in), optional :: range
    integer :: i
    do i = 1, size(p)
       call make_random_momentum (p(i), range)
    end do
  end subroutine make_random_momentum_vector

  subroutine make_random_vector (v, range)
    type(vector), intent(inout) :: v
    real(kind=default), intent(in), optional :: range
    call make_random_complex (v%t, range)
    call make_random_complex_vector (v%x, range)
  end subroutine make_random_vector

  subroutine make_random_vector_vector (v, range)
    type(vector), dimension(:), intent(inout) :: v
    real(kind=default), intent(in), optional :: range
    integer :: i
    do i = 1, size(v)
       call make_random_vector (v(i), range)
    end do
  end subroutine make_random_vector_vector

  subroutine make_random_spinor (psi, range)
    type(spinor), intent(inout) :: psi
    real(kind=default), intent(in), optional :: range
    call make_random_complex_vector (psi%a, range)
  end subroutine make_random_spinor

  subroutine make_random_spinor_vector (psi, range)
    type(spinor), dimension(:), intent(inout) :: psi
    real(kind=default), intent(in), optional :: range
    integer :: i
    do i = 1, size(psi)
       call make_random_spinor (psi(i), range)
    end do
  end subroutine make_random_spinor_vector

  subroutine make_random_conjspinor (psibar, range)
    type(conjspinor), intent(inout) :: psibar
    real(kind=default), intent(in), optional :: range
    call make_random_complex_vector (psibar%a, range)
  end subroutine make_random_conjspinor

  subroutine make_random_conjspinor_vector (psibar, range)
    type(conjspinor), dimension(:), intent(inout) :: psibar
    real(kind=default), intent(in), optional :: range
    integer :: i
    do i = 1, size(psibar)
       call make_random_conjspinor (psibar(i), range)
    end do
  end subroutine make_random_conjspinor_vector

  subroutine make_random_tensor (t, range)
    type(tensor), intent(inout) :: t
    real(kind=default), intent(in), optional :: range
    call make_random_complex_array (t%t, range)
  end subroutine make_random_tensor

  subroutine make_random_tensor_vector (t, range)
    type(tensor), dimension(:), intent(inout) :: t
    real(kind=default), intent(in), optional :: range
    integer :: i
    do i = 1, size(t)
       call make_random_tensor (t(i), range)
    end do
  end subroutine make_random_tensor_vector

  subroutine make_random_tensor2odd (t, range)
    type(tensor2odd), intent(inout) :: t
    real(kind=default), intent(in), optional :: range
    call make_random_complex_vector (t%e, range)
    call make_random_complex_vector (t%b, range)
  end subroutine make_random_tensor2odd

  subroutine make_random_tensor2odd_vector (t, range)
    type(tensor2odd), dimension(:), intent(inout) :: t
    real(kind=default), intent(in), optional :: range
    integer :: i
    do i = 1, size(t)
       call make_random_tensor2odd (t(i), range)
    end do
  end subroutine make_random_tensor2odd_vector

  subroutine make_random_bispinor (chi, range)
    type(bispinor), intent(inout) :: chi
    real(kind=default), intent(in), optional :: range
    call make_random_complex_vector (chi%a, range)
  end subroutine make_random_bispinor

  subroutine make_random_bispinor_vector (chi, range)
    type(bispinor), dimension(:), intent(inout) :: chi
    real(kind=default), intent(in), optional :: range
    integer :: i
    do i = 1, size(chi)
       call make_random_bispinor (chi(i), range)
    end do
  end subroutine make_random_bispinor_vector

end module keystones_tools
