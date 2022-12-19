! test_openmp.f90 --
! driver.f90 -- O'Mega self test driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!      Christian Speckner <cnspeckn@googlemail.com>
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

program test_openmp

  use omp_lib
  use parameters_QCD, init_parameters_qcd => init_parameters

  use kinds
  use constants
  use amplitude_openmp
  use omega_testtools

  implicit none

  integer, parameter :: NOUT = 3
  integer, parameter :: NCALLS = 2000
  real(kind=default), parameter :: ROOTS = 1000

  real(kind=default) :: elapsed
  real(kind=default), dimension(:), allocatable :: amp2
  integer :: max_threads, num_threads
  logical :: passed

  call init_parameters_qcd ()

  passed = .true.

  max_threads = omp_get_max_threads ()
  write (unit = *, fmt = "(1X, 'max. threads: ', I3)")  max_threads
  allocate (amp2(0:max_threads))

  call omp_set_dynamic (.true.)
  call evaluate (NCALLS, ROOTS, 1, 1, amp2(0), elapsed)
  write (unit = *, fmt = "(1X, A, F8.4, A, F8.4, A, E19.10)") &
       '      dynamic: elapsed ', elapsed, ' seconds, elapsed * #threads: ', &
       elapsed * max_threads, ' seconds, amp2 = ', amp2(0)

  call omp_set_dynamic (.false.)
  do num_threads = 1, max_threads
     call omp_set_num_threads (num_threads)
     call evaluate (NCALLS, ROOTS, 1, 1, amp2(num_threads), elapsed)
     write (unit = *, fmt = "(1X, A, I2, A, F8.4, A, F8.4, A, E19.10)") &
          '#threads = ', num_threads, ', elapsed ', elapsed, &
          ' seconds, elapsed * #threads: ', elapsed*num_threads, &
          ' seconds, amp2 = ', amp2(num_threads)
  end do

  call expect (minval(amp2), maxval(amp2), "test_openmp", passed, &
       quiet=.false., threshold = 0.60_default)

  if (passed) then
     stop 0
  else
     stop 1
  end if

  contains

    subroutine evaluate (n, roots, flv, hel, amp2, elapsed)
      integer, intent(in) :: n
      real(kind=default), intent(in) :: roots
      integer, intent(in) :: flv, hel
      real(kind=default), intent(out) :: amp2, elapsed

      real(kind=default) :: wtime_start, wtime
      real(kind=default) :: sum_amp2
      real(kind=default) :: p(0:3,2+NOUT)
      integer :: i, size
      integer, dimension(:), allocatable :: seed

      call random_seed (size)
      allocate (seed(size))
      seed = 42
      call random_seed (put = seed)
      deallocate (seed)

      call beams (roots, 0.0_default, 0.0_default, p(:,1), p(:,2))
      sum_amp2 = 0

      wtime_start = omp_get_wtime ()
      do i = 1, n
         call massless_isotropic_decay (roots, p(:,3:))
         call new_event (p)
         sum_amp2 = sum_amp2 + color_sum (1, 1)
      end do
      elapsed = omp_get_wtime () - wtime_start

      amp2 = sum_amp2 / n

    end subroutine evaluate

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
      real(kind=default), dimension(4) :: ran
      real(kind=default) :: c, s, f, qabs, x, r, z
      integer :: k
      ! Generate isotropic null vectors
      do k = 1, size (p, dim = 2)
         call random_number (ran)
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

end program test_openmp

