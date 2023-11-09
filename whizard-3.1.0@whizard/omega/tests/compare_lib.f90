! compare_lib.f90 -- compare two O'Mega versions
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

module compare_lib
  ! use ieee_arithmetic
  use kinds
  use constants
  use tao_random_numbers
  use omega95
  use omega_interface
  use omega_testtools
  implicit none
  private
  public :: check
  public :: omega_flavor_states, omega_squared_matrix_element
  public :: massless_isotropic_decay, rambo, beams, dot, rambo_check
contains

  elemental function ieee_is_nan (x) result (yorn)
    logical :: yorn
    real (kind=default), intent(in) :: x
    yorn = (x /= x)
  end function ieee_is_nan

  subroutine check (v1, v2, roots, threshold, n, &
                    failures, attempts, seed, abs_threshold, &
                    ignore_phase, flip_sign)
    type(omega_procedures), intent(in) :: v1, v2
    real(kind=default), intent(in) :: roots, threshold
    integer, intent(in) :: n
    integer, intent(out) :: failures, attempts
    integer, intent(in), optional :: seed
    real(kind=default), intent(in), optional :: abs_threshold
    logical, intent(in), optional :: ignore_phase, flip_sign
    logical :: modulus_only
    logical :: match, passed
    integer :: n_out, n_flv, n_hel, n_col
    integer :: i, i_flv, i_hel, i_col
    real(kind=default), dimension(:,:), allocatable :: p
    complex(kind=default) :: a1, a2
    real(kind=default) :: asq1, asq2, s_asq1, s_asq2, relative_sign
    character(len=80) :: msg
    modulus_only = .false.
    if (present (ignore_phase)) then
       modulus_only = ignore_phase
    end if
    relative_sign = 1
    if (present (flip_sign)) then
       if (flip_sign) then
          relative_sign = -1
       end if
    end if
    failures = 0
    attempts = 0
    a1 = 0
    a2 = 0
    asq1 = 0
    asq2 = 0
    s_asq1 = 0
    s_asq2 = 0
    call quantum_numbers (v1, v2, n_out, n_flv, n_hel, n_col, match)
    if (.not.match) then
       failures = 1
       return
    end if
    if (n_out <= 0) then
       print *, "no outgoing particles"
       failures = 1
       return
    end if
    if (n_flv <= 0) then
       print *, "no allowed flavor combinations"
       failures = 1
       return
    end if
    if (n_hel <= 0) then
       print *, "no allowed helicity combinations"
       failures = 1
       return
    end if
    if (n_col <= 0) then
       print *, "no allowed color flows"
       failures = 1
       return
    end if
    if (present (seed)) then
       call tao_random_seed (seed)
    end if
    call v1%reset_helicity_selection (-1.0_default, -1)
    call v2%reset_helicity_selection (-1.0_default, -1)
    allocate (p(0:3,2+n_out))
    call beams (ROOTS, 0.0_default, 0.0_default, p(:,1), p(:,2))
    do i = 1, N
       if (n_out > 1) then
          call massless_isotropic_decay (ROOTS, p(:,3:))
       end if
       if (n_out == 1) then
          p(:,3) = p(:,1) + p(:,2)
       end if
       call v1%new_event (p)
       call v2%new_event (p)
       do i_flv = 1, n_flv
          do i_hel = 1, n_hel
             attempts = attempts + 1
             passed = .true.
             do i_col = 1, n_col
                a1 = v1%get_amplitude (i_flv, i_hel, i_col)
                a2 = v2%get_amplitude (i_flv, i_hel, i_col)*relative_sign
                if (ieee_is_nan (real (a1)) .or. ieee_is_nan (aimag (a1))) then
                   write (*, "(1X,'evt=',I5,', flv=',I3,', col=',I3,': ', A)") &
                        i, i_flv, i_col, "v1 amplitude NaN"
                end if
                if (ieee_is_nan (real (a2)) .or. ieee_is_nan (aimag (a2))) then
                   write (*, "(1X,'evt=',I5,', flv=',I3,', col=',I3,': ', A)") &
                        i, i_flv, i_col, "v2 amplitude NaN"
                end if
                write (msg, "(1X,'evt=',I5,', flv=',I3,', col=',I3,', hel=',I3)") &
                     i, i_flv, i_col, i_hel
                if (modulus_only) then
                   call expect (abs (a1), abs (a2), trim(msg), passed, &
                                quiet=.true., threshold=threshold, &
                                abs_threshold=abs_threshold)
                else
                   call expect (a1, a2, trim(msg), passed, &
                                quiet=.true., threshold=threshold, &
                                abs_threshold=abs_threshold)
                end if
             end do
             write (msg, "(1X,'evt=',I5,', flv=',I3,', hel=',I3)") &
                  i, i_flv, i_hel
             asq1 = v1%color_sum (i_flv, i_hel)
             s_asq1 = s_asq1 + asq1
             asq2 = v2%color_sum (i_flv, i_hel)
             s_asq2 = s_asq2 + asq2
             call expect (asq1, asq2, trim(msg), passed, &
                          quiet=.true., threshold=threshold, &
                          abs_threshold=abs_threshold)
             if (.not.passed) then
                failures = failures + 1
             end if
          end do
       end do
    end do
    print *, 'Summed results: '
    print *, 's_asq1, s_asq2 =    ', s_asq1, s_asq2
    deallocate (p)
  end subroutine check

  subroutine quantum_numbers (v1, v2, n_out, n_flv, n_hel, n_col, match)
    type(omega_procedures), intent(in) :: v1, v2
    integer, intent(out) :: n_out, n_flv, n_hel, n_col
    logical, intent(out) :: match
    integer, dimension(:,:), allocatable :: &
         v1_flavor_states, v2_flavor_states, &
         v1_spin_states, v2_spin_states
    integer, dimension(:,:,:), allocatable :: &
         v1_color_flows, v2_color_flows
    logical, dimension(:,:), allocatable :: &
         v1_ghost_flags, v2_ghost_flags
    type(omega_color_factor), dimension(:), allocatable :: &
         v1_color_factors, v2_color_factors
    integer :: n_in, n_prt, n_cix, n_cfs
    n_in = v1%number_particles_in ()
    n_out = v1%number_particles_out ()
    n_prt = n_in + n_out
    n_flv = v1%number_flavor_states ()
    n_hel = v1%number_spin_states ()
    n_cix = v1%number_color_indices ()
    n_col = v1%number_color_flows ()
    n_cfs = v1%number_color_factors ()
    match = .true.
    if (v2%number_particles_in () .ne. n_in) then
       print *, "number_particles_in don't match!"
       match = .false.
    end if
    if (v2%number_particles_out () .ne. n_out) then
       print *, "number_particles_out don't match!"
       match = .false.
    end if
    if (v2%number_flavor_states () .ne. n_flv) then
       print *, "number_flavor_states don't match!"
       match = .false.
    end if
    if (v2%number_spin_states () .ne. n_hel) then
       print *, "number_spin_states don't match!"
       match = .false.
    end if
    if (v2%number_color_indices () .ne. n_cix) then
       print *, "number_color_indices don't match!"
       match = .false.
    end if
    if (v2%number_color_flows () .ne. n_col) then
       print *, "number_color_flows don't match!"
       match = .false.
    end if
    ! We save only the symmetric part in the OVM
    !if (v2%number_color_factors () .ne. n_cfs) then
       !print *, "number_color_factors don't match!"
       !match = .false.
    !end if
    if (match) then
       allocate (v1_flavor_states(n_prt,n_flv), v2_flavor_states(n_prt,n_flv))
       allocate (v1_spin_states(n_prt,n_hel), v2_spin_states(n_prt,n_hel))
       allocate (v1_color_flows(n_cix,n_prt,n_col), &
                 v2_color_flows(n_cix,n_prt,n_col))
       allocate (v1_ghost_flags(n_prt,n_col), v2_ghost_flags(n_prt,n_col))
       !allocate (v1_color_factors(n_cfs), v2_color_factors(n_cfs))
       call v1%flavor_states (v1_flavor_states)
       call v2%flavor_states (v2_flavor_states)
       call v1%spin_states (v1_spin_states)
       call v2%spin_states (v2_spin_states)
       call v1%color_flows (v1_color_flows, v1_ghost_flags)
       call v2%color_flows (v2_color_flows, v2_ghost_flags)
       !call v1%color_factors (v1_color_factors)
       !call v2%color_factors (v2_color_factors)
       if (any (v1_flavor_states .ne. v2_flavor_states)) then
          print *, "flavor states don't match!"
          print *, "CAVEAT: this might be due to simple reordering!"
          match = .false.
       end if
       if (any (v1_spin_states .ne. v2_spin_states)) then
          print *, "spin states don't match!"
          print *, "CAVEAT: this might be due to simple reordering!"
          match = .false.
       end if
       if (any (v1_color_flows .ne. v2_color_flows)) then
          print *, "color flows don't match!"
          print *, "CAVEAT: this might be due to simple reordering!"
          match = .false.
       end if
       if (any (v1_ghost_flags .neqv. v2_ghost_flags)) then
          print *, "ghost flags don't match!"
          print *, "CAVEAT: this might be due to simple reordering!"
          match = .false.
       end if
       !if (any (.not. color_factors_equal (v1_color_factors, &
                                           !v2_color_factors))) then
          !print *, "color_factors don't match!"
          !print *, "CAVEAT: this might be due to simple reordering!"
          !match = .false.
       !end if
       deallocate (v1_flavor_states, v2_flavor_states)
       deallocate (v1_spin_states, v2_spin_states)
       deallocate (v1_color_flows, v2_color_flows)
       deallocate (v1_ghost_flags, v2_ghost_flags)
       !deallocate (v1_color_factors, v2_color_factors)
    end if
  end subroutine quantum_numbers

  elemental function color_factors_equal (cf1, cf2) result (eq)
    logical :: eq
    type(omega_color_factor), intent(in) :: cf1, cf2
    eq = (cf1%i1 .eq. cf2%i1) .and. (cf1%i2 .eq. cf2%i2) .and. (cf1%factor .eq. cf2%factor)
  end function color_factors_equal

  subroutine omega_flavor_states (proc, flavors)
    type(omega_procedures) :: proc
    integer, dimension(:,:), allocatable, intent(inout) :: flavors
    integer :: n_in, n_out, n_prt, n_flv
    n_in = proc%number_particles_in ()
    n_out = proc%number_particles_out ()
    n_prt = n_in + n_out
    n_flv = proc%number_flavor_states ()
    if (allocated (flavors)) then
       if (any (size (flavors) /= (/ n_prt, n_flv /))) then
          deallocate (flavors)
          allocate (flavors (n_prt, n_flv))
       end if
    else
       allocate (flavors (n_prt, n_flv))
    end if
    call proc%flavor_states (flavors)
  end subroutine omega_flavor_states

  subroutine omega_squared_matrix_element (proc, p, asq, error)
    type(omega_procedures) :: proc
    real(kind=default), dimension(0:,:), intent(in) :: p
    real(kind=default), intent(out) :: asq
    logical, intent(out) :: error
    real(kind=default) :: asq_sum
    integer :: i_hel
    call proc%new_event (p)
    error = .false.
    if (proc%number_flavor_states () /= 1) then
       print *, "ambiguous flavor in omega amplitude"
       error = .true.
       return
    end if
    asq_sum = 0
    do i_hel = 1, proc%number_spin_states ()
       asq_sum = asq_sum + proc%color_sum (1, i_hel)
    end do
    asq = asq_sum / 4
  end subroutine omega_squared_matrix_element

  pure function dot (p, q) result (pq)
    real(kind=default), dimension(0:), intent(in) :: p, q
    real(kind=default) :: pq
    pq = p(0)*q(0) - dot_product (p(1:), q(1:))
  end function dot

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

  !------------------------------------------------------
  ! RAMBO, R. Kleiss, W.J. Stirling, S.D. Ellis.
  ! Comp. Phys. Commun. 40 (1986) 359
  !------------------------------------------------------
  subroutine rambo (roots, m, p, weight, unweighted)
    implicit none
    real(kind=default), intent(in) :: roots
    real(kind=default), dimension(:), intent(in) :: m
    real(kind=default), dimension(0:,:), intent(out) :: p
    real(kind=default), intent(out) :: weight
    logical, intent(in) :: unweighted

    real(kind=default), dimension(0:3,size(m)):: q
    real(kind=default), dimension(0:3) :: sum_q
    real(kind=default) :: mass_sum_q

    real(kind=default), dimension(4) :: random
    real(kind=default) :: random_weight
    real(kind=double), dimension(4) :: random_double
    real(kind=double) :: random_weight_double

    real(kind=default), dimension(size(m)):: m2, e, v, p2

    real(kind=default) :: a,accu,bq,costh,phi,f0,g,g0,pm2
    real(kind=default) :: sinth,sm2,w,wt2,wt3,wtm,wtmax,x,x2,xmax,sum_m
    real(kind=default) :: b(3)

    integer :: i, iter, k, num_massive

    real(kind=default), dimension(:), allocatable, save :: z
    real(kind=default), save :: twopi, log_pi_over_2
    real(kind=default), parameter :: ACC = 1d-14

    integer, parameter :: MAX_ITERATIONS = 6
    integer, save :: underflows = 0, overflows = 0, excessive_weights = 0

    if (size(p,dim=2) /= size(m)) then
       print *, 'rambo: mismatch of array dimensions of M and P'
       stop
    end if

    ! initialize the factorials for the phase space weight
    if (allocated(z)) then
       if (size(z) < size(m)) then
          deallocate (z)
       end if
    end if
    if (.not.allocated(z)) then
       allocate (z(size(m)))
       ! z(1) = ???
       twopi = 8 * atan (1.0_default)
       log_pi_over_2 = log (twopi / 4)
       z(2) = log_pi_over_2
       do k = 3, size(z)
          z(k) = z(k-1) + log_pi_over_2 - 2 * log (real (k-2, kind=default))
       end do
       do k = 3, size(z)
          z(k) = z(k) - log (real (k-1, kind=default))
       end do
    end if

    ! check whether total energy suffices and count nonzero masses
    num_massive = count (m /= 0)
    sum_m = sum (abs (m))
    if (sum_m > roots) then
       print *, ' RAMBO FAILS: TOTAL MASS =', sum_m, &
                ' IS NOT', ' SMALLER THAN TOTAL ENERGY =', roots
       stop
    end if

    ! generate N massless momenta
    generate: do

       do i = 1, size(m)
          call tao_random_number (random_double)
          random = random_double
          costh = 2 * random(1) - 1
          sinth = sqrt (1 - costh*costh)
          phi = twopi * random(2)
          q(0,i) = -log (random(3)*random(4))
          q(3,i) = q(0,i) * costh
          q(2,i) = q(0,i) * sinth * cos (phi)
          q(1,i) = q(0,i) * sinth * sin (phi)
       end do

       ! compute the parameters of the conformal transformation
       sum_q = sum (q, dim=2)
       mass_sum_q = sqrt (dot (sum_q, sum_q))
       b = - sum_q(1:3) / mass_sum_q
       g = sum_q(0) / mass_sum_q
       a = 1 / (1 + g)
       x = roots / mass_sum_q

       ! TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
       do i = 1, size(m)
          bq = b(1) * q(1,i) + b(2) * q(2,i) + b(3) * q(3,i)
          p(1:3,i) = x * (q(1:3,i) + b * (q(0,i) + a*bq))
          p(0,i) = x * (g*q(0,i) + bq)
       end do

       ! for unweighted massless momenta, we're done
       weight = 1
       if (num_massive == 0 .and. unweighted) then
          exit generate
       end if

       ! CALCULATE WEIGHT AND POSSIBLE WARNINGS
       weight = log_pi_over_2
       if (size(m) /= 2) then
          weight = (2*size(m) - 4) * log (roots) + z(size(m))
       end if
       if (weight < - 180) then
          if (underflows <= 5) then
             call rambo_flow (weight, 'under')
          end if
          underflows = underflows + 1
       end if
       if (weight > 174) then
          if (overflows <= 5) then
             call rambo_flow (weight, 'over')
          end if
          overflows = overflows + 1
       end if

       ! return FOR WEIGHTED MASSLESS MOMENTA
       if (num_massive /= 0) then

          ! MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
          xmax = sqrt (1 - (sum_m/roots)**2)
          m2 = m**2
          p2 = p(0,:)**2

          x = xmax
          accu = roots * ACC

          iter = 0
          solve: do
             f0 = - roots
             g0 = 0
             x2 = x*x
             do i = 1, size(m)
                e(i) = sqrt (m2(i) + x2 * p2(i))
                f0 = f0 + e(i)
                g0 = g0 + p2(i) / e(i)
             end do
             if (abs (f0) > accu) then
                iter = iter + 1
                if (iter <= MAX_ITERATIONS) then
                   x = x - f0 / (x*g0)
                   cycle solve
                else
                   print *, ' RAMBO WARNS:', MAX_ITERATIONS, &
                        ' ITERATIONS DID NOT GIVE THE', &
                        ' DESIRED ACCURACY =', ACC
                end if
             end if
             exit solve
          end do solve

          v = x * p(0,:)
          P(1:3,:) = x * P(1:3,:)
          P(0,:) = e

          ! CALCULATE THE MASS-EFFECT WEIGHT FACTOR
          wt2 = product (v / e)
          wt3 = sum (v**2 / e)
          wtm = (2*size(m) - 3) * log (x) + LOG (wt2 / wt3 * roots)

          if (unweighted) then

             ! UNWEIGHTED MASSIVE MOMENTA REQUIRED: ESTIMATE MAXIMUM WEIGHT
             weight = exp (wtm)
             if (num_massive <= 1) then
                ! ONE MASSIVE PARTICLE
                wtmax = xmax**(4*size(m) - 6)
             elseif (num_massive > 2) then
                ! MORE THAN TWO MASSIVE PARTICLES: AN ESTIMATE ONLY
                wtmax = xmax**(2*size(m) - 5 + num_massive)
             else
                ! TWO MASSIVE PARTICLES
                sm2 = sum (m2)
                ! this was wrong (always 0) in thr orignal)
                pm2 = product (m2, mask = (m2 /= 0))
                wtmax = ((1 - sm2 / (roots**2))**2 &
                            - 4*pm2 / roots**4)**(size(m) - 1.5_default)
             end if

             ! DETERMINE WHETHER OR NOT TO ACCEPT THIS EVENT
             w = weight / wtmax
             if (w > 1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the purpose of our tests, we can suppress this warning
!               print *, ' RAMBO WARNS: ESTIMATE FOR MAXIMUM WEIGHT =', &
!                        wtmax, '    EXCEEDED BY A FACTOR ', w
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                excessive_weights = excessive_weights + 1
             end if
             call tao_random_number (random_weight_double)
             random_weight = random_weight_double
             if (w < random_weight) then
                cycle generate
             end if
             weight = 1

          else

             ! return FOR  WEIGHTED MASSIVE MOMENTA
             weight = weight + wtm
             if (weight < -180) then
                if (underflows <= 5) then
                   call rambo_flow (weight, 'under')
                end if
                underflows = underflows + 1
             end if
             if (weight > 174) then
                if (overflows <= 5) then
                   call rambo_flow (weight, 'over')
                end if
                overflows = overflows + 1
             end if
             weight = exp (weight)
          end if

       else
          weight = exp (weight)
       end if

       exit generate
    end do generate

  end subroutine rambo

  subroutine rambo_check (roots, m, p, quiet)
    real(kind=default), intent(in) :: roots
    real(kind=default), dimension(:), intent(in) :: m
    real(kind=default), dimension(0:,:), intent(in) :: p
    logical, intent (in) :: quiet
    real(kind=default), dimension(0:3) :: sum_p
    integer :: mu, i
    logical :: passed
    real(kind=default), parameter :: &
         THRESHOLD_MOMENTUM = 0.80, &
         THRESHOLD_MASS = 0.45
    passed = .true.
    sum_p = sum (p, dim=2)
    call expect (sum_p(0), roots, 'energy momentum', &
                 passed, threshold=THRESHOLD_MOMENTUM, quiet=quiet)
    do mu = 1, 3
       call expect_zero (sum_p(mu), roots, 'spatial momentum', &
                         passed, threshold=THRESHOLD_MOMENTUM, quiet=quiet)
    end do
    do i = 1, size(m)
       call expect (dot (p(:,i), p(:,i)), m(i)**2, 'mass shell', &
                         passed, threshold=THRESHOLD_MASS, quiet=quiet)
    end do
    if (.not.passed .and. .not.quiet) then
       do i = 1, size (m)
          print *, 'M(', i, ') = ', sqrt (abs (dot (p(:,i), p(:,i)))), &
                   'vs. ', m(i)
       end do
       do mu = 0, 3
          print *, 'sum p(', mu, ',:) = ', sum_p(mu)
       end do
    end if
  end subroutine rambo_check

  subroutine rambo_flow (w, f)
    implicit none
    real(kind=default), intent(in) :: w
    character(len=*), intent(in) :: f
    print *, ' RAMBO WARNS: WEIGHT = EXP(', w,') MAY ', f
  end subroutine rambo_flow
    
end module compare_lib
