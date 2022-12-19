! fermi_lib.f90 -- check On Shell Ward Identities in O'Mega
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

module fermi_lib
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
contains

  elemental function ieee_is_nan (x) result (yorn)
    logical :: yorn
    real (kind=default), intent(in) :: x
    yorn = (x /= x)
  end function ieee_is_nan

  subroutine check (process, i, j, eps, roots, threshold, abs_threshold, &
       n, failures, attempts, seed)
    type(omega_procedures), intent(in) :: process
    integer, intent(in) :: i, j, eps
    real(kind=default), intent(in) :: roots, threshold, abs_threshold
    integer, intent(in) :: n
    integer, intent(out) :: failures, attempts
    integer, intent(in), optional :: seed
    logical :: match, passed
    integer :: n_in, n_out, n_prt, n_flv, n_hel, n_col, n_cix
    integer :: k, k_flv, k_hel, k_col, k_xhl, k_xcl
    ! integer :: k_prt
    integer, dimension(:,:,:), allocatable :: color_states
    logical, dimension(:,:), allocatable :: ghost_flags
    integer, dimension(:,:), allocatable :: spin_states
    integer, dimension(:), allocatable :: spin_map, color_map
    real(kind=default), dimension(:,:), allocatable :: p
    complex(kind=default), dimension(:,:,:), allocatable :: a1, a2
    complex(kind=default) :: aa1, aa2
    character(len=80) :: msg
    failures = 0
    attempts = 0
    if (present (seed)) then
       call tao_random_seed (seed)
    end if
    n_in = process%number_particles_in ()
    n_out = process%number_particles_out ()
    n_prt = n_in + n_out
    n_flv = process%number_flavor_states ()
    n_hel = process%number_spin_states ()
    n_col = process%number_color_flows ()
    n_cix = process%number_color_indices ()
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
    if (max (i, j) > n_prt .or. min (i, j) < 1) then
       print *, "invalid #particles, i or j!"
       stop 2
    end if
    call process%reset_helicity_selection (-1.0_default, -1)
    allocate (p(0:3,n_prt))
    allocate (a1(n_flv,n_hel,n_col), a2(n_flv,n_hel,n_col))
    allocate (spin_states(n_prt,n_hel))
    allocate (spin_map(n_hel))
    allocate (color_states(n_cix,n_prt,n_col))
    allocate (ghost_flags(n_prt,n_col))
    allocate (color_map(n_col))
    call process%spin_states(spin_states)
    call spin_exchange (spin_map, spin_states, i, j)
    call process%color_flows (color_states,ghost_flags)
    ! call color_exchange (color_map, color_states, i, j)
    if (n_in == 2) then
       call beams (ROOTS, 0.0_default, 0.0_default, p(:,1), p(:,2))
    else if (n_in == 1) then
       p(0,1) = ROOTS
       p(1:3,1) = 0
    else
       stop
    end if
    do k = 1, N
       attempts = attempts + 1
       passed = .true.
       call massless_isotropic_decay (ROOTS, p(:,3:))
       call process%new_event (p)
       do k_flv = 1, n_flv
          do k_col = 1, n_col
             do k_hel = 1, n_hel
                a1(k_flv,k_hel,k_col) &
                     = process%get_amplitude (k_flv, k_hel, k_col)
             end do
          end do
       end do
       ! print *, p
       call exchange_momentum (p, i, j)
       ! print *, p
       call process%new_event (p)
       do k_flv = 1, n_flv
          do k_col = 1, n_col
             do k_hel = 1, n_hel
                a2(k_flv,k_hel,k_col) &
                     = process%get_amplitude (k_flv, k_hel, k_col)
             end do
          end do
       end do
       do k_flv = 1, n_flv
          do k_col = 1, n_col
             k_xcl = k_col
             ! k_xcl = color_map(k_col)
             ! print *, 'c (', k_col, ') = ', color_states(:,:,k_col)
             ! print *, 'cx(', k_xcl, ') = ', color_states(:,:,k_xcl)
             do k_hel = 1, n_hel
                k_xhl = spin_map(k_hel)
                ! print *, 'h  = ', spin_states(:,k_hel)
                ! print *, 'hx = ', spin_states(:,k_xhl)
                aa1 = a1(k_flv,k_hel,k_col)
                aa2 = a2(k_flv,k_xhl,k_xcl)
                if (ieee_is_nan (real (aa1)) .or. ieee_is_nan (aimag (aa1))) then
                   write (*, "(1X,'evt=',I5,', flv=',I3, &
                               &', hel=',I3,', col=',I3,': ', A)") &
                        k, k_flv, k_hel, k_col, "v1 amplitude NaN"
                end if
                if (ieee_is_nan (real (aa2)) .or. ieee_is_nan (aimag (aa2))) then
                   write (*, "(1X,'evt=',I5,', flv=',I3, &
                               &', hel=',I3,', col=',I3,': ', A)") &
                        k, k_flv, k_xhl, k_xcl, "v2 amplitude NaN"
                end if
                write (msg, "(1X,'evt=',I5,', flv=',I3, &
                              &', col=',I3,',',I3,', hel=',I3,',',I3)") &
                     k, k_flv, k_col, k_xcl, k_hel, k_xhl
                call expect (aa1, eps * aa2, trim(msg), passed, &
                                quiet=.true., threshold=threshold, &
                                abs_threshold=abs_threshold)
                if (.not.passed) then
                   failures = failures + 1
                end if
             end do
          end do
       end do
    end do
    deallocate (p)
    deallocate (a1, a2)
    deallocate (spin_states)
  end subroutine check

  subroutine exchange_momentum (p, i, j)
    real(kind=default), dimension(0:,:), intent(inout) :: p
    integer, intent(in) :: i, j
    real(kind=default), dimension(0:ubound(p,1)) :: tmp
    tmp = p(:,j)
    p(:,j) = p(:,i)
    p(:,i) = tmp
  end subroutine exchange_momentum
  
  subroutine exchange_spins (h, i, j)
    integer, dimension(:,:), intent(inout) :: h
    integer, intent(in) :: i, j
    integer, dimension(size(h,2)) :: tmp
    tmp = h(j,:)
    h(j,:) = h(i,:)
    h(i,:) = tmp
  end subroutine exchange_spins
  
  subroutine spin_exchange (map, h, i, j)
    integer, dimension(:), intent(inout) :: map
    integer, dimension(:,:), intent(in) :: h
    integer, intent(in) :: i, j
    integer, dimension(size(h,1),size(h,2)) :: hx
    integer :: k, l
    hx = h
    call exchange_spins (hx, i, j)
    do l = 1, size (h, 2)
       find: do k = 1, size (h, 2)
          if (all (hx(:,k) == h(:,l))) then
             map(l) = k
             exit find
          endif
       end do find
    end do
  end subroutine spin_exchange

  subroutine exchange_colors (c, i, j)
    integer, dimension(:,:,:), intent(inout) :: c
    integer, intent(in) :: i, j
    integer, dimension(size(c,1),size(c,3)) :: tmp
    tmp = c(:,j,:)
    c(:,j,:) = c(:,i,:)
    c(:,i,:) = tmp
  end subroutine exchange_colors

  subroutine color_exchange (map, c, i, j)
    integer, dimension(:), intent(inout) :: map
    integer, dimension(:,:,:), intent(in) :: c
    integer, intent(in) :: i, j
    integer, dimension(size(c,1),size(c,2),size(c,3)) :: cx
    integer :: k, l
    cx = c
    call exchange_colors (cx, i, j)
    do l = 1, size (c, 3)
       print *, 'c (', l, ') = ', c(:,:,l)
    end do
    do l = 1, size (c, 3)
       print *, 'cx(', l, ') = ', cx(:,:,l)
    end do
    do l = 1, size (c, 3)
       map(l) = l
       find: do k = 1, size (c, 3)
          ! this does NOT work!!!
          ! color flow indices can be renamed
          ! for equivalent flows ...
          if (all (cx(:,:,k) == c(:,:,l))) then
             print *, 'map: ', l, ' -> ', k
             map(l) = k
             exit find
          endif
       end do find
    end do
  end subroutine color_exchange
  
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

end module fermi_lib
