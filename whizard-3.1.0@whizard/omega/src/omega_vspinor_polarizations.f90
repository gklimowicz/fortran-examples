!  omegalib.nw --
!
!  Copyright (C) 1999-2022 by
!      Wolfgang Kilian <kilian@physik.uni-siegen.de>
!      Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!      Juergen Reuter <juergen.reuter@desy.de>
!      with contributions from                                                                                                                                    
!      Fabian Bach <fabian.bach@t-online.de>                                                                                                                 
!      Bijan Chokoufe Nejad <bijan.chokoufe@desy.de>                                                                                                              
!      Christian Speckner <cnspeckn@googlemail.com>     
!
!  WHIZARD is free software; you can redistribute it and/or modify it
!  under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2, or (at your option)
!  any later version.
!
!  WHIZARD is distributed in the hope that it will be useful, but
!  WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module omega_vspinor_polarizations
  use kinds
  use constants
  use omega_vectors
  use omega_bispinors
  use omega_bispinor_couplings
  use omega_vectorspinors
  implicit none
  public :: ueps, veps
  private :: eps
  private :: outer_product
  integer, parameter, public :: omega_vspinor_pols_2010_01_A = 0
contains
  pure function eps (mass, k, s) result (e)
    type(vector) :: e
    real(kind=default), intent(in) :: mass
    type(momentum), intent(in) :: k
    integer, intent(in) :: s
    real(kind=default) :: kabs, kabs2, sqrt2, m
    real(kind=default) :: cos_phi, sin_phi, cos_th, sin_th
    complex(kind=default) :: epiphi, emiphi
    sqrt2 = sqrt (2.0_default)
    kabs2 = dot_product (k%x, k%x)
    m = abs(mass)
    if (kabs2 > 0) then
       kabs = sqrt (kabs2)
       if ((k%x(1) == 0) .and. (k%x(2) == 0)) then
          cos_phi = 1
          sin_phi = 0
       else
          cos_phi = k%x(1) / sqrt(k%x(1)**2 + k%x(2)**2)
          sin_phi = k%x(2) / sqrt(k%x(1)**2 + k%x(2)**2)
       end if
       cos_th = k%x(3) / kabs
       sin_th = sqrt(1 - cos_th**2)
       epiphi = cos_phi + (0,1) * sin_phi
       emiphi = cos_phi - (0,1) * sin_phi
       e%t = 0
       e%x = 0
       select case (s)
       case (1)
          e%x(1) = epiphi * (-cos_th * cos_phi + (0,1) * sin_phi) / sqrt2
          e%x(2) = epiphi * (-cos_th * sin_phi - (0,1) * cos_phi) / sqrt2
          e%x(3) = epiphi * ( sin_th / sqrt2)
       case (-1)
          e%x(1) = emiphi * ( cos_th * cos_phi + (0,1) * sin_phi) / sqrt2
          e%x(2) = emiphi * ( cos_th * sin_phi - (0,1) * cos_phi) / sqrt2
          e%x(3) = emiphi * (-sin_th / sqrt2)
       case (0)
          if (m > 0) then
             e%t = kabs / m
             e%x = k%t / (m*kabs) * k%x
          end if
       case (4)
          if (m > 0) then
             e = (1 / m) * k
          else
             e = (1 / k%t) * k
          end if
       end select
    else   !!! for particles in their rest frame defined to be
           !!! polarized along the 3-direction
       e%t = 0
       e%x = 0
       select case (s)
       case (1)
          e%x(1) = cmplx ( - 1,   0, kind=default) / sqrt2
          e%x(2) = cmplx (   0,   1, kind=default) / sqrt2
       case (-1)
          e%x(1) = cmplx (   1,   0, kind=default) / sqrt2
          e%x(2) = cmplx (   0,   1, kind=default) / sqrt2
       case (0)
          if (m > 0) then
             e%x(3) = 1
          end if
       case (4)
          if (m > 0) then
             e = (1 / m) * k
          else
             e = (1 / k%t) * k
          end if
       end select
    end if
  end function eps
  pure function ueps (m, k, s) result (t)
    type(vectorspinor) :: t
    real(kind=default), intent(in) :: m
    type(momentum), intent(in) :: k
    integer, intent(in) :: s
    integer :: i
    type(vector) :: ep, e0, em
    type(bispinor) :: up, um
    do i = 1, 4
      t%psi(i)%a = 0
    end do
    select case (s)
    case (2)
       ep = eps (m, k, 1)
       up = u (m, k, 1)
       t = outer_product (ep, up)
    case (1)
       ep = eps (m, k, 1)
       e0 = eps (m, k, 0)
       up = u (m, k, 1)
       um = u (m, k, -1)
       t = (1 / sqrt (3.0_default)) * (outer_product (ep, um) &
            + sqrt (2.0_default) * outer_product (e0, up))
    case (-1)
       e0 = eps (m, k, 0)
       em = eps (m, k, -1)
       up = u (m, k, 1)
       um = u (m, k, -1)
       t = (1 / sqrt (3.0_default)) * (sqrt (2.0_default) * &
            outer_product (e0, um) + outer_product (em, up))
    case (-2)
       em = eps (m, k, -1)
       um = u (m, k, -1)
       t = outer_product (em, um)
    end select
  end function ueps
  pure function veps (m, k, s) result (t)
    type(vectorspinor) :: t
    real(kind=default), intent(in) :: m
    type(momentum), intent(in) :: k
    integer, intent(in) :: s
    integer :: i
    type(vector) :: ep, e0, em
    type(bispinor) :: vp, vm
    do i = 1, 4
      t%psi(i)%a = 0
    end do
    select case (s)
    case (2)
       ep = conjg(eps (m, k, 1))
       vp = v (m, k, 1)
       t = outer_product (ep, vp)
    case (1)
       ep = conjg(eps (m, k, 1))
       e0 = conjg(eps (m, k, 0))
       vp = v (m, k, 1)
       vm = v (m, k, -1)
       t = (1 / sqrt (3.0_default)) * (outer_product (ep, vm) &
            + sqrt (2.0_default) * outer_product (e0, vp))
    case (-1)
       e0 = conjg(eps (m, k,  0))
       em = conjg(eps (m, k, -1))
       vp = v (m, k, 1)
       vm = v (m, k, -1)
       t = (1 / sqrt (3.0_default)) * (sqrt (2.0_default) &
            * outer_product (e0, vm) + outer_product (em, vp))
    case (-2)
       em = conjg(eps (m, k, -1))
       vm = v (m, k, -1)
       t = outer_product (em, vm)
    end select
  end function veps
  pure function outer_product (ve, sp) result (vs)
    type(vectorspinor) :: vs
    type(vector), intent(in) :: ve
    type(bispinor), intent(in) :: sp
    integer :: i
    vs%psi(1)%a(1:4) = ve%t * sp%a(1:4)
    do i = 1, 3
       vs%psi((i+1))%a(1:4) = ve%x(i) * sp%a(1:4)
    end do
  end function outer_product
end module omega_vspinor_polarizations
