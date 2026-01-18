! kinematics.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module kinematics
  use kinds
  use constants
  use products, only: dot
  use specfun, only: gamma
  implicit none
  private
  public :: boost_velocity
  private :: boost_one_velocity, boost_many_velocity
  public :: boost_momentum
  private :: boost_one_momentum, boost_many_momentum
  public :: lambda
  public :: two_to_three
  private :: two_to_three_massive, two_to_three_massless
  public :: one_to_two
  private :: one_to_two_massive, one_to_two_massless
  public :: polar_to_cartesian, on_shell
  public :: massless_isotropic_decay
  public :: phase_space_volume
  interface boost_velocity
     module procedure boost_one_velocity, boost_many_velocity
  end interface
  interface boost_momentum
     module procedure boost_one_momentum, boost_many_momentum
  end interface
  interface two_to_three
     module procedure two_to_three_massive, two_to_three_massless
  end interface
  interface one_to_two
     module procedure one_to_two_massive, one_to_two_massless
  end interface
  type, public :: LIPS3
     real(kind=default), dimension(3,0:3) :: p
     real(kind=default) :: jacobian
  end type LIPS3
contains
  pure function boost_one_velocity (p, beta) result (p_prime)
    real(kind=default), dimension(0:), intent(in) :: p
    real(kind=default), dimension(1:), intent(in) :: beta
    real(kind=default), dimension(0:3) :: p_prime
    real(kind=default), dimension(1:3) :: b
    real(kind=default) :: gamma, b_dot_p
    gamma = 1.0 / sqrt (1.0 - dot_product (beta, beta))
    b = gamma * beta
    b_dot_p = dot_product (b, p(1:3))
    p_prime(0) = gamma * p(0) - b_dot_p
    p_prime(1:3) = p(1:3) + (b_dot_p / (1.0 + gamma) - p(0)) * b
  end function boost_one_velocity
  pure function boost_many_velocity (p, beta) result (p_prime)
    real(kind=default), dimension(:,0:), intent(in) :: p
    real(kind=default), dimension(1:), intent(in) :: beta
    real(kind=default), dimension(size(p,dim=1),0:3) :: p_prime
    integer :: i
    do i = 1, size (p, dim=1)
       p_prime(i,:) = boost_one_velocity (p(i,:), beta)
    end do
  end function boost_many_velocity
  pure function boost_one_momentum (p, q) result (p_prime)
    real(kind=default), dimension(0:), intent(in) :: p, q
    real(kind=default), dimension(0:3) :: p_prime
    p_prime = boost_velocity (p, q(1:3) / abs (q(0)))
  end function boost_one_momentum
  pure function boost_many_momentum (p, q) result (p_prime)
    real(kind=default), dimension(:,0:), intent(in) :: p
    real(kind=default), dimension(0:), intent(in) :: q
    real(kind=default), dimension(size(p,dim=1),0:3) :: p_prime
    p_prime = boost_many_velocity (p, q(1:3) / abs (q(0)))
  end function boost_many_momentum
  pure function lambda (a, b, c) result (lam)
    real(kind=default), intent(in) :: a, b, c
    real(kind=default) :: lam
    lam = a**2 + b**2 + c**2 - 2*(a*b + b*c + c*a)
  end function lambda
  pure function two_to_three_massive &
       (s, t1, s2, phi, cos_theta3, phi3, ma, mb, m1, m2, m3) result (p)
    real(kind=default), intent(in) :: &
         s, t1, s2, phi, cos_theta3, phi3, ma, mb, m1, m2, m3
    type(LIPS3) :: p
    real(kind=default), dimension(0:3) :: p23
    real(kind=default) :: Ea, pa_abs, E1, p1_abs, p3_abs, cos_theta
    pa_abs = sqrt (lambda (s, ma**2, mb**2) / (4 * s))
    Ea = sqrt (ma**2 + pa_abs**2)
    p1_abs = sqrt (lambda (s, m1**2, s2) / (4 * s))
    E1 = sqrt (m1**2 + p1_abs**2)
    p3_abs = sqrt (lambda (s2, m2**2, m3**2) / (4 * s2))
    p%jacobian = &
         1.0 / (2*PI)**5 * (p3_abs / pa_abs) / (32 * sqrt (s * s2))
    cos_theta = (t1 - ma**2 - m1**2 + 2*Ea*E1) / (2*pa_abs*p1_abs)
    p%p(1,1:3) = polar_to_cartesian (p1_abs, cos_theta, phi)
    p%p(1,0) = on_shell (p%p(1,:), m1)
    p23(1:3) = - p%p(1,1:3)
    p23(0) = on_shell (p23, sqrt (s2))
    p%p(3:2:-1,:) = one_to_two (p23, cos_theta3, phi3, m3, m2)
  end function two_to_three_massive
  pure function two_to_three_massless (s, t1, s2, phi, cos_theta3, phi3) &
       result (p)
    real(kind=default), intent(in) :: s, t1, s2, phi, cos_theta3, phi3
    type(LIPS3) :: p
    real(kind=default), dimension(0:3) :: p23
    real(kind=default) :: pa_abs, p1_abs, p3_abs, cos_theta
    pa_abs = sqrt (s) / 2
    p1_abs = (s - s2) / (2 * sqrt (s))
    p3_abs = sqrt (s2) / 2
    p%jacobian = 1.0 / ((2*PI)**5 * 32 * s)
    cos_theta = 1 + t1 / (2*pa_abs*p1_abs)
    p%p(1,0) = p1_abs
    p%p(1,1:3) = polar_to_cartesian (p1_abs, cos_theta, phi)
    p23(1:3) = - p%p(1,1:3)
    p23(0) = on_shell (p23, sqrt (s2))
    p%p(3:2:-1,:) = one_to_two (p23, cos_theta3, phi3)
  end function two_to_three_massless
  pure function one_to_two_massive (p12, cos_theta, phi, m1, m2) result (p)
    real(kind=default), dimension(0:), intent(in) :: p12
    real(kind=default), intent(in) :: cos_theta, phi, m1, m2
    real(kind=default), dimension(2,0:3) :: p
    real(kind=default) :: s, p1_abs
    s = dot (p12, p12)
    p1_abs = sqrt (lambda (s, m1**2, m2**2) / (4 * s))
    p(1,1:3) = polar_to_cartesian (p1_abs, cos_theta, phi)
    p(2,1:3) = - p(1,1:3)
    p(1,0) = on_shell (p(1,:), m1)
    p(2,0) = on_shell (p(2,:), m2)
    p = boost_momentum (p, - p12)
  end function one_to_two_massive
  pure function one_to_two_massless (p12, cos_theta, phi) result (p)
    real(kind=default), dimension(0:), intent(in) :: p12
    real(kind=default), intent(in) :: cos_theta, phi
    real(kind=default), dimension(2,0:3) :: p
    real(kind=default) :: p1_abs
    p1_abs = sqrt (dot (p12, p12)) / 2
    p(1,0) = p1_abs
    p(1,1:3) = polar_to_cartesian (p1_abs, cos_theta, phi)
    p(2,0) = p1_abs
    p(2,1:3) = - p(1,1:3)
    p = boost_momentum (p, - p12)
  end function one_to_two_massless
  pure function polar_to_cartesian (v_abs, cos_theta, phi) result (v)
    real(kind=default), intent(in) :: v_abs, cos_theta, phi
    real(kind=default), dimension(3) :: v
    real(kind=default) :: sin_phi, cos_phi, sin_theta
    sin_theta = sqrt (1.0 - cos_theta**2)
    cos_phi = cos (phi)
    sin_phi = sin (phi)
    v = (/ sin_theta * cos_phi, sin_theta * sin_phi, cos_theta /) * v_abs
  end function polar_to_cartesian
  pure function on_shell (p, m) result (E)
    real(kind=default), dimension(0:), intent(in) :: p
    real(kind=default), intent(in) :: m
    real(kind=default) :: E
    E = sqrt (m**2 + dot_product (p(1:3), p(1:3)))
  end function on_shell
  pure function massless_isotropic_decay (roots, ran) result (p)
    real (kind=default), intent(in) :: roots
    real (kind=default), dimension(:,:), intent(in) :: ran
    real (kind=default), dimension(size(ran,dim=1),0:3) :: p
    real (kind=default), dimension(size(ran,dim=1),0:3) :: q
    real (kind=default), dimension(0:3) :: qsum
    real (kind=default) :: cos_theta, sin_theta, phi, qabs, x, r, z
    integer :: k
    do k = 1, size (p, dim = 1)
       q(k,0) = - log (ran(k,1) * ran(k,2))
       cos_theta = 2 * ran(k,3) - 1
       sin_theta = sqrt (1 - cos_theta**2)
       phi = 2 * PI * ran(k,4)
       q(k,1) = q(k,0) * sin_theta * cos (phi)
       q(k,2) = q(k,0) * sin_theta * sin (phi)  
       q(k,3) = q(k,0) * cos_theta
    enddo
    qsum = sum (q, dim = 1)
    qabs = sqrt (dot (qsum, qsum))
    x = roots / qabs
    do k = 1, size (p, dim = 1)
       r = dot (q(k,:), qsum) / qabs
       z = (q(k,0) + r) / (qsum(0) + qabs)
       p(k,1:3) = x * (q(k,1:3) - qsum(1:3) * z)
       p(k,0) = x * r
    enddo
  end function massless_isotropic_decay
  pure function phase_space_volume (n, roots) result (volume)
    integer, intent(in) :: n
    real (kind=default), intent(in) :: roots
    real (kind=default) :: volume
    real (kind=default) :: nd
    nd = n
    volume = (nd - 1) / (8*PI * (gamma (nd))**2) * (roots / (4*PI))**(2*n-4)
  end function phase_space_volume
end module kinematics
module phase_space
  use kinds
  use constants
  use kinematics !NODEP!
  use tao_random_numbers
  implicit none
  private
  public :: random_LIPS3
  private :: random_LIPS3_unit, random_LIPS3_unit_massless
  private :: LIPS3_unit_to_s2_t1_angles, LIPS3_unit_to_s2_t1_angles_m0
  interface random_LIPS3
     module procedure random_LIPS3_unit, random_LIPS3_unit_massless
  end interface
  type, public :: LIPS3_unit
     real(kind=default), dimension(5) :: x
     real(kind=default) :: s
     real(kind=default), dimension(2) :: mass_in
     real(kind=default), dimension(3) :: mass_out
     real(kind=default) :: jacobian
  end type LIPS3_unit
  type, public :: LIPS3_unit_massless
     real(kind=default), dimension(5) :: x
     real(kind=default) :: s
     real(kind=default) :: jacobian
  end type LIPS3_unit_massless
  type, public :: LIPS3_s2_t1_angles
     real(kind=default) :: s2, t1, phi, cos_theta3, phi3
     real(kind=default) :: s
     real(kind=default), dimension(2) :: mass_in
     real(kind=default), dimension(3) :: mass_out
     real(kind=default) :: jacobian
  end type LIPS3_s2_t1_angles
  type, public :: LIPS3_s2_t1_angles_massless
     real(kind=default) :: s2, t1, phi, cos_theta3, phi3
     real(kind=default) :: s
     real(kind=default) :: jacobian
  end type LIPS3_s2_t1_angles_massless
  type, public :: LIPS3_momenta
     real(kind=default), dimension(0:3,3) :: p
     real(kind=default) :: s
     real(kind=default), dimension(2) :: mass_in
     real(kind=default), dimension(3) :: mass_out
     real(kind=default) :: jacobian
  end type LIPS3_momenta
  type, public :: LIPS3_momenta_massless
     real(kind=default), dimension(0:3,3) :: p
     real(kind=default) :: s
     real(kind=default) :: jacobian
  end type LIPS3_momenta_massless
contains
  pure subroutine random_LIPS3_unit (rng, lips)
    type(tao_random_state), intent(inout) :: rng
    type(LIPS3_unit), intent(inout) :: lips
    call tao_random_number (rng, lips%x)
    lips%jacobian = 1
  end subroutine random_LIPS3_unit
  pure subroutine random_LIPS3_unit_massless (rng, lips)
    type(tao_random_state), intent(inout) :: rng
    type(LIPS3_unit_massless), intent(inout) :: lips
    call tao_random_number (rng, lips%x)
    lips%jacobian = 1
  end subroutine random_LIPS3_unit_massless
  pure subroutine LIPS3_unit_to_s2_t1_angles (s2_t1_angles, unit)
    type(LIPS3_s2_t1_angles), intent(out) :: s2_t1_angles
    type(LIPS3_unit), intent(in) :: unit
  end subroutine  LIPS3_unit_to_s2_t1_angles
  pure subroutine LIPS3_unit_to_s2_t1_angles_m0 (s2_t1_angles, unit)
    type(LIPS3_s2_t1_angles_massless), intent(out) :: s2_t1_angles
    type(LIPS3_unit_massless), intent(in) :: unit
  end subroutine  LIPS3_unit_to_s2_t1_angles_m0
end module phase_space
