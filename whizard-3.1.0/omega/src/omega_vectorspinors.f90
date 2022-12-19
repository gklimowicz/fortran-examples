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
module omega_vectorspinors
  use kinds
  use constants
  use omega_bispinors
  use omega_vectors
  implicit none
  private
  public :: operator (*), operator (+), operator (-)
  public :: abs, set_zero
  type, public :: vectorspinor
     ! private (omegalib needs access, but DON'T TOUCH IT!)
     type(bispinor), dimension(4) :: psi
  end type vectorspinor
  interface operator (*)
    module procedure vspinor_product
  end interface
  private :: vspinor_product
  interface set_zero
    module procedure set_zero_vectorspinor
  end interface
  private :: set_zero_vectorspinor
  interface operator (*)
     module procedure integer_vectorspinor, vectorspinor_integer, &
            real_vectorspinor, double_vectorspinor, &
            complex_vectorspinor, dcomplex_vectorspinor, &
            vectorspinor_real, vectorspinor_double, &
            vectorspinor_complex, vectorspinor_dcomplex, &
            momentum_vectorspinor, vectorspinor_momentum
  end interface
  private :: integer_vectorspinor, vectorspinor_integer, real_vectorspinor, &
       double_vectorspinor, complex_vectorspinor, dcomplex_vectorspinor, &
       vectorspinor_real, vectorspinor_double, vectorspinor_complex, &
       vectorspinor_dcomplex
  interface operator (+)
     module procedure plus_vectorspinor
  end interface
  private :: plus_vectorspinor
  interface operator (-)
     module procedure neg_vectorspinor
  end interface
  private :: neg_vectorspinor
  interface operator (+)
     module procedure add_vectorspinor
  end interface
  private :: add_vectorspinor
  interface operator (-)
     module procedure sub_vectorspinor
  end interface
  private :: sub_vectorspinor
  interface abs
     module procedure abs_vectorspinor
  end interface
  private :: abs_vectorspinor
  integer, parameter, public :: omega_vectorspinors_2010_01_A = 0
contains
  pure function vspinor_product (psil, psir) result (psilpsir)
    complex(kind=default) :: psilpsir
    type(vectorspinor), intent(in) :: psil, psir
    psilpsir = psil%psi(1) * psir%psi(1) &
             - psil%psi(2) * psir%psi(2) &
             - psil%psi(3) * psir%psi(3) &
             - psil%psi(4) * psir%psi(4)
  end function vspinor_product
  elemental subroutine set_zero_vectorspinor (x)
    type(vectorspinor), intent(out) :: x
    call set_zero (x%psi)
  end subroutine set_zero_vectorspinor
  pure function integer_vectorspinor (x, y) result (xy)
    type(vectorspinor) :: xy
    integer, intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
      xy%psi(k) = x * y%psi(k)
    end do
  end function integer_vectorspinor
  pure function real_vectorspinor (x, y) result (xy)
    type(vectorspinor) :: xy
    real(kind=single), intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
    xy%psi(k) = x * y%psi(k)
    end do
  end function real_vectorspinor
  pure function double_vectorspinor (x, y) result (xy)
    type(vectorspinor) :: xy
    real(kind=default), intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
    xy%psi(k) = x * y%psi(k)
    end do
  end function double_vectorspinor
  pure function complex_vectorspinor (x, y) result (xy)
    type(vectorspinor) :: xy
    complex(kind=single), intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
    xy%psi(k) = x * y%psi(k)
    end do
  end function complex_vectorspinor
  pure function dcomplex_vectorspinor (x, y) result (xy)
    type(vectorspinor) :: xy
    complex(kind=default), intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
    xy%psi(k) = x * y%psi(k)
    end do
  end function dcomplex_vectorspinor
  pure function vectorspinor_integer (y, x) result (xy)
    type(vectorspinor) :: xy
    integer, intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
    xy%psi(k) = y%psi(k) * x
    end do
  end function vectorspinor_integer
  pure function vectorspinor_real (y, x) result (xy)
    type(vectorspinor) :: xy
    real(kind=single), intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
    xy%psi(k) = y%psi(k) * x
    end do
  end function vectorspinor_real
  pure function vectorspinor_double (y, x) result (xy)
    type(vectorspinor) :: xy
    real(kind=default), intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
    xy%psi(k) = y%psi(k) * x
    end do
  end function vectorspinor_double
  pure function vectorspinor_complex (y, x) result (xy)
    type(vectorspinor) :: xy
    complex(kind=single), intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
    xy%psi(k) = y%psi(k) * x
    end do
  end function vectorspinor_complex
  pure function vectorspinor_dcomplex (y, x) result (xy)
    type(vectorspinor) :: xy
    complex(kind=default), intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
    xy%psi(k) = y%psi(k) * x
    end do
  end function vectorspinor_dcomplex
  pure function momentum_vectorspinor (y, x) result (xy)
    type(bispinor) :: xy
    type(momentum), intent(in) :: y
    type(vectorspinor), intent(in) :: x
    integer :: k
    do k = 1,4
    xy%a(k) = y%t    * x%psi(1)%a(k) - y%x(1) * x%psi(2)%a(k) - &
            y%x(2) * x%psi(3)%a(k) - y%x(3) * x%psi(4)%a(k)
    end do
  end function momentum_vectorspinor
  pure function vectorspinor_momentum (y, x) result (xy)
    type(bispinor) :: xy
    type(momentum), intent(in) :: x
    type(vectorspinor), intent(in) :: y
    integer :: k
    do k = 1,4
    xy%a(k) = x%t    * y%psi(1)%a(k) - x%x(1) * y%psi(2)%a(k) - &
            x%x(2) * y%psi(3)%a(k) - x%x(3) * y%psi(4)%a(k)
    end do
  end function vectorspinor_momentum
  pure function plus_vectorspinor (x) result (plus_x)
    type(vectorspinor) :: plus_x
    type(vectorspinor), intent(in) :: x
    integer :: k
    do k = 1,4
    plus_x%psi(k) = + x%psi(k)
    end do
  end function plus_vectorspinor
  pure function neg_vectorspinor (x) result (neg_x)
    type(vectorspinor) :: neg_x
    type(vectorspinor), intent(in) :: x
    integer :: k
    do k = 1,4
    neg_x%psi(k) = - x%psi(k)
    end do
  end function neg_vectorspinor
  pure function add_vectorspinor (x, y) result (xy)
    type(vectorspinor) :: xy
    type(vectorspinor), intent(in) :: x, y
    integer :: k
    do k = 1,4
    xy%psi(k) = x%psi(k) + y%psi(k)
    end do
  end function add_vectorspinor
  pure function sub_vectorspinor (x, y) result (xy)
    type(vectorspinor) :: xy
    type(vectorspinor), intent(in) :: x, y
    integer :: k
    do k = 1,4
    xy%psi(k) = x%psi(k) - y%psi(k)
    end do
  end function sub_vectorspinor
  pure function abs_vectorspinor (psi) result (x)
    real(kind=default) :: x
    type(vectorspinor), intent(in) :: psi
    x = sqrt (real (dot_product (psi%psi(1)%a, psi%psi(1)%a) &
            - dot_product (psi%psi(2)%a, psi%psi(2)%a)   &
            - dot_product (psi%psi(3)%a, psi%psi(3)%a)   &
            - dot_product (psi%psi(4)%a, psi%psi(4)%a)))
  end function abs_vectorspinor
end module omega_vectorspinors
