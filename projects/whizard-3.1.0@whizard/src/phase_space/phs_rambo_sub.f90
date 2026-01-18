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

submodule (phs_rambo) phs_rambo_s

  use io_units
  use constants
  use numeric_utils
  use format_defs, only: FMT_19
  use permutations, only: factorial
  use diagnostics
  use physics_defs

  implicit none

  integer, parameter :: BISECT_MAX_ITERATIONS = 1000
  real(default), parameter :: BISECT_MIN_PRECISION = tiny_10

contains

  module subroutine phs_rambo_config_final (object)
    class(phs_rambo_config_t), intent(inout) :: object
  end subroutine phs_rambo_config_final

  module subroutine phs_rambo_config_write (object, unit, include_id)
    class(phs_rambo_config_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: include_id
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Partonic, flat phase-space configuration (RAMBO):"
    call object%base_write (unit)
  end subroutine phs_rambo_config_write

  module subroutine phs_rambo_config_configure (phs_config, sqrts, &
       sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, &
       ignore_mismatch, nlo_type, subdir)
    class(phs_rambo_config_t), intent(inout) :: phs_config
    real(default), intent(in) :: sqrts
    logical, intent(in), optional :: sqrts_fixed
    logical, intent(in), optional :: lab_is_cm
    logical, intent(in), optional :: azimuthal_dependence
    logical, intent(in), optional :: rebuild
    logical, intent(in), optional :: ignore_mismatch
    integer, intent(in), optional :: nlo_type
    type(string_t), intent(in), optional :: subdir
    if (.not. present (nlo_type)) &
      phs_config%nlo_type = BORN
    if (phs_config%n_out < 2) then
       call msg_fatal ("RAMBO phase space requires n_out >= 2")
    end if
    phs_config%n_channel = 1
    phs_config%n_par = 3 * phs_config%n_out - 4
    phs_config%sqrts = sqrts
    if (present (sqrts_fixed))  phs_config%sqrts_fixed = sqrts_fixed
    if (present (lab_is_cm))  phs_config%lab_is_cm = lab_is_cm
    if (allocated (phs_config%channel))  deallocate (phs_config%channel)
    allocate (phs_config%channel (1))
    call phs_config%compute_md5sum ()
  end subroutine phs_rambo_config_configure

  module subroutine phs_rambo_config_startup_message (phs_config, unit)
    class(phs_rambo_config_t), intent(in) :: phs_config
    integer, intent(in), optional :: unit
    call phs_config%base_startup_message (unit)
    write (msg_buffer, "(A,2(1x,I0,1x,A))") &
         "Phase space: flat (RAMBO)"
    call msg_message (unit = unit)
  end subroutine phs_rambo_config_startup_message

  module subroutine phs_rambo_write (object, unit, verbose)
    class(phs_rambo_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    u = given_output_unit (unit)
    call object%base_write (u)
    write (u, "(1X,A)") "Intermediate masses (massless):"
    write (u, "(3X,999(" // FMT_19 // "))") object%k
    write (u, "(1X,A)") "Intermediate masses (massive):"
    write (u, "(3X,999(" // FMT_19 // "))") object%m
  end subroutine phs_rambo_write

  module subroutine phs_rambo_final (object)
    class(phs_rambo_t), intent(inout) :: object
  end subroutine phs_rambo_final

  module subroutine phs_rambo_init (phs, phs_config)
    class(phs_rambo_t), intent(out) :: phs
    class(phs_config_t), intent(in), target :: phs_config
    call phs%base_init (phs_config)
    associate (n => phs%config%n_out)
      select case (n)
      case (1)
         if (sum (phs%m_out) > phs%m_in (1)) then
            print *, "m_in = ", phs%m_in
            print *, "m_out = ", phs%m_out
            call msg_fatal &
                 ("[phs_rambo_init] Decay is kinematically forbidden.")
         end if
      end select
      allocate (phs%k(n), source = 0._default)
      allocate (phs%m(n), source = 0._default)
      phs%volume = 1. / (twopi)**(3 * n) &
           * (pi / 2.)**(n - 1) / (factorial(n - 1) * factorial(n - 2))
    end associate
  end subroutine phs_rambo_init

  module subroutine phs_rambo_evaluate_selected_channel (phs, c_in, r_in)
    class(phs_rambo_t), intent(inout) :: phs
    integer, intent(in) :: c_in
    real(default), intent(in), dimension(:) :: r_in
    type(vector4_t), dimension(2) :: p_rest, p_boosted
    type(vector4_t) :: q
    real(default), dimension(2) :: r_angle
    integer :: i
    if (.not. phs%p_defined) return
    call phs%select_channel (c_in)
    phs%r(:,c_in) = r_in
    associate (n => phs%config%n_out, m => phs%m)
      call phs%generate_intermediates (r_in(:n - 2))
      q = sum (phs%p)
      do i = 2, n
         r_angle(1) = r_in(n - 5 + 2 * i)
         r_angle(2) = r_in(n - 4 + 2 * i)
         call phs%decay_intermediate (i, r_angle, p_rest)
         p_boosted = boost(q, m(i - 1)) * p_rest
         q = p_boosted(1)
         phs%q(i - 1)  = p_boosted(2)
      end do
      phs%q(n) = q
    end associate
    phs%q_defined = .true.
    phs%r_defined = .true.
  end subroutine phs_rambo_evaluate_selected_channel

  module subroutine phs_rambo_evaluate_other_channels (phs, c_in)
    class(phs_rambo_t), intent(inout) :: phs
    integer, intent(in) :: c_in
  end subroutine phs_rambo_evaluate_other_channels

  module subroutine phs_rambo_decay_intermediate (phs, i, r_angle, p)
    class(phs_rambo_t), intent(in) :: phs
    integer, intent(in) :: i
    real(default), dimension(2), intent(in) :: r_angle
    type(vector4_t), dimension(2), intent(out) :: p
    real(default) :: k_abs, cos_theta, phi
    type(vector3_t):: k
    real(default), dimension(2) :: E
    cos_theta = 2. * r_angle(1) - 1.
    phi = twopi * r_angle(2)
    if (phi > pi) phi = phi - twopi
    k_abs = sqrt (lambda (phs%m(i - 1)**2, phs%m(i)**2, phs%m_out(i - 1)**2)) &
         / (2. * phs%m(i - 1))
    k = k_abs * [cos(phi) * sqrt(1. - cos_theta**2), &
         sin(phi) * sqrt(1. - cos_theta**2), cos_theta]
    E(1) = sqrt (phs%m(i)**2 + k_abs**2)
    E(2) = sqrt (phs%m_out(i - 1)**2 + k_abs**2)
    p(1) = vector4_moving (E(1), -k)
    p(2) = vector4_moving (E(2), k)
  end subroutine phs_rambo_decay_intermediate

  module subroutine phs_rambo_generate_intermediates (phs, r)
    class(phs_rambo_t), intent(inout) :: phs
    real(default), dimension(:), intent(in) :: r
    integer :: i, j
    associate (n => phs%config%n_out, k => phs%k, m => phs%m)
      m(1) = invariant_mass (sum (phs%p))
      m(n) = phs%m_out (n)
      call calculate_k (r)
      do i = 2, n - 1
         m(i) = k(i) + sum (phs%m_out (i:n))
      end do
      ! Massless volume times reweighting for massive volume
      phs%f(1) = k(1)**(2 * n - 4) &
           * 8. * rho(m(n - 1), phs%m_out(n), phs%m_out(n - 1))
      do i = 2, n - 1
         phs%f(1) = phs%f(1) * &
              rho(m(i - 1), m(i), phs%m_out(i - 1)) / &
              rho(k(i - 1), k(i), 0._default) * &
              M(i) / K(i)
      end do
    end associate
  contains
    subroutine calculate_k (r)
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), allocatable :: u
      integer :: i
      associate (n => phs%config%n_out, k => phs%k, m => phs%m)
        k = 0
        k(1) = m(1) - sum(phs%m_out(1:n))
        allocate (u(2:n - 1), source=0._default)
        call solve_for_u (r, u)
        do i = 2, n - 1
           k(i) = sqrt (u(i) * k(i - 1)**2)
        end do
      end associate
    end subroutine calculate_k

    subroutine solve_for_u (r, u)
      real(default), dimension(phs%config%n_out - 2), intent(in) :: r
      real(default), dimension(2:phs%config%n_out - 1), intent(out) :: u
      integer :: i, j
      real(default) :: f, f_mid, xl, xr, xmid
      associate (n => phs%config%n_out)
        do i = 2, n - 1
           xl = 0
           xr = 1
           if (r(i - 1) == 1 .or. r(i - 1) == 0) then
              u(i) = r(i - 1)
           else
              do j = 1, BISECT_MAX_ITERATIONS
                 xmid = (xl + xr) / 2.
                 f = f_rambo (xl, n - i) - r(i - 1)
                 f_mid = f_rambo (xmid, n - i) - r(i - 1)
                 if (f * f_mid > 0) then
                    xl = xmid
                 else
                    xr = xmid
                 end if
                 if (abs(xl - xr) < BISECT_MIN_PRECISION) exit
              end do
              u(i) = xmid
           end if
        end do
      end associate
    end subroutine solve_for_u

    real(default) function f_rambo(u, n)
      real(default), intent(in) :: u
      integer, intent(in) :: n
      f_rambo = (n + 1) * u**n - n * u**(n + 1)
    end function f_rambo

    real(default) function rho (M1, M2, m)
      real(default), intent(in) :: M1, M2, m
      real(default) :: MP, MM
      rho = sqrt ((M1**2 - (M2 + m)**2) * (M1**2 - (M2 - m)**2))
      ! MP = (M1 - (M2 + m)) * (M1 + (M2 + m))
      ! MM = (M1 - (M2 - m)) * (M1 + (M2 - m))
      ! rho = sqrt (MP) * sqrt (MM)
      rho = rho / (8._default * M1**2)
    end function rho

  end subroutine phs_rambo_generate_intermediates

  module subroutine phs_rambo_invert_intermediates (phs)
    class(phs_rambo_t), intent(inout) :: phs
    real(default) :: u
    integer :: i
    associate (n => phs%config%n_out, k => phs%k, m => phs%m)
      k = m
      do i = 1, n - 1
         k(i) = k(i) - sum (phs%m_out(i:n))
      end do
      do i = 2, n - 1
         u = (k(i) / k(i - 1))**2
         phs%r(i - 1, 1) = (n + 1 - i) * u**(n - i) &
              - (n - i) * u**(n + 1 - i)
      end do
    end associate
  end subroutine phs_rambo_invert_intermediates

  module subroutine phs_rambo_inverse (phs)
    class(phs_rambo_t), intent(inout) :: phs
    type(vector4_t), dimension(:), allocatable :: q
    type(vector4_t) :: p
    type(lorentz_transformation_t) :: L
    real(default) :: phi, cos_theta
    integer :: i
    if (.not. (phs%p_defined .and. phs%q_defined)) return
    call phs%select_channel ()
    associate (n => phs%config%n_out, m => phs%m)
      allocate(q(n))
      m(1) = invariant_mass (sum (phs%p))
      q(1) = vector4_at_rest (m(1))
      q(n) = phs%q(n)
      do i = 2, n - 1
         q(i) = q(i) + sum (phs%q(i:n))
         m(i) = invariant_mass (q(i))
      end do
      call phs%invert_intermediates ()
      do i = 2, n
         L = inverse (boost (q(i - 1), m(i - 1)))
         p = L * phs%q(i - 1)
         phi = azimuthal_angle (p); cos_theta = polar_angle_ct (p)
         phs%r(n - 5 + 2 * i, 1) = (cos_theta + 1.) / 2.
         phs%r(n - 4 + 2 * i, 1) = phi / twopi
      end do
    end associate
       phs%r_defined = .true.
  end subroutine phs_rambo_inverse


end submodule phs_rambo_s

