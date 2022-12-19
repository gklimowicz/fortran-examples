! benchmark_lib.f90 -- benchmark two O'Mega versions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 1999-2014 by
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

module benchmark_lib
  ! use ieee_arithmetic
  use kinds
  use constants
  use tao_random_numbers
  use omega95
  use omega_interface
  use omega_testtools
  implicit none
  private
  public :: benchmark, beams, massless_isotropic_decay

contains

  subroutine benchmark (time1, time2, v1, v2, roots, n, seed)
    type(omega_procedures), intent(in) :: v1, v2
    real(double), intent(out) :: time1, time2
    real(kind=default), intent(in) :: roots
    integer, intent(in) :: n
    integer, intent(in), optional :: seed
    real(kind=double) :: wtime_start, wtime
    integer :: n_out
    integer :: i
    real(kind=default), dimension(:,:), allocatable :: p
    n_out = v1%number_particles_out ()
    call v1%reset_helicity_selection (-1.0_default, -1)
    call v2%reset_helicity_selection (-1.0_default, -1)
    allocate (p(0:3,2+n_out))
    call beams (ROOTS, 0.0_default, 0.0_default, p(:,1), p(:,2))

    if (present (seed)) then
       call tao_random_seed (seed)
    end if
    call cpu_time (wtime_start)
    do i = 1, n
       call massless_isotropic_decay (roots, p(:,3:))
       call v1%new_event (p)
    end do
    call cpu_time (wtime)
    time1 = 1000 * (wtime - wtime_start) / n

    ! Use the same seed for both to compare on equal terms
    if (present (seed)) then
       call tao_random_seed (seed)
    end if
    call cpu_time (wtime_start)
    do i = 1, n
       call massless_isotropic_decay (roots, p(:,3:))
       call v2%new_event (p)
    end do
    call cpu_time (wtime)
    time2 = 1000 * (wtime - wtime_start) / n

  end subroutine benchmark

  pure function dot (p, q) result (pq)
    real(kind=default), dimension(0:), intent(in) :: p, q
    real(kind=default) :: pq
    pq = p(0)*q(0) - dot_product (p(1:), q(1:))
  end function dot

  pure function mass2 (p) result (m2)
    real(kind=default), dimension(0:), intent(in) :: p
    real(kind=default) :: m2
    m2 = p(0)*p(0) - p(1)*p(1) - p(2)*p(2) - p(3)*p(3)
  end function mass2

  pure subroutine beams (roots, m1, m2, p1, p2)
    real(kind=default), intent(in) :: roots, m1, m2
    real(kind=default), dimension(0:), intent(out) :: p1, p2
    real(kind=default) :: m12, m22
    m12 = m1**2
    m22 = m2**2
    p1(0) = (roots**2 + m12 - m22) / (2*roots)
    p1(1:2) = 0
    p1(3) = sqrt (p1(0)**2 - m12)
    p2(0) = roots - p1(0)
    p2(1:3) = - p1(1:3)
  end subroutine beams

  ! The massless RAMBO algorithm
  subroutine massless_isotropic_decay (roots, p)
    real(kind=default), intent(in) :: roots
    real(kind=default), dimension(0:,:), intent(out) :: p
    real(kind=default), dimension(0:3,size(p,dim=2)) :: q
    real(kind=default), dimension(0:3) :: qsum
    real(kind=double), dimension(4) :: ran_double
    real(kind=default), dimension(4) :: ran
    real(kind=default) :: c, s, f, qabs, x, r, z
    integer :: k
    ! Generate isotropic null vectors
    do k = 1, size (p, dim = 2)
       ! if default is not double or single, we can't use
       ! tao_random_number directly ...
       call tao_random_number (ran_double)
       ran = ran_double
       ! generate a x*exp(-x) distribution for q(0,k)
       q(0,k)= -log(ran(1)*ran(2))
       c = 2*ran(3)-1
       f = 2*PI*ran(4)
       s = sqrt(1-c*c)
       q(2,k) = q(0,k)*s*sin(f)
       q(3,k) = q(0,k)*s*cos(f)
       q(1,k) = q(0,k)*c
    enddo
    ! Boost and rescale the vectors
    qsum = sum (q, dim = 2)
    qabs = sqrt (dot (qsum, qsum))
    x = roots/qabs
    do k = 1, size (p, dim = 2)
       r = dot (q(0:,k), qsum) / qabs
       z = (q(0,k)+r)/(qsum(0)+qabs)
       p(1:3,k) = x*(q(1:3,k)-qsum(1:3)*z)
       p(0,k) = x*r
    enddo
  end subroutine massless_isotropic_decay

end module benchmark_lib
