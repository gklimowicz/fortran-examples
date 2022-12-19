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
module omega_spinors
  use kinds
  use constants
  implicit none
  private
  public :: operator (*), operator (+), operator (-)
  public :: abs, set_zero
  
  type, public :: conjspinor
     ! private (omegalib needs access, but DON'T TOUCH IT!)
     complex(kind=default), dimension(4) :: a
  end type conjspinor
  type, public :: spinor
     ! private (omegalib needs access, but DON'T TOUCH IT!)
     complex(kind=default), dimension(4) :: a
  end type spinor
  interface operator (*)
     module procedure conjspinor_spinor
  end interface
  private :: conjspinor_spinor
  interface set_zero
    module procedure set_zero_spinor, set_zero_conjspinor
  end interface
  private :: set_zero_spinor, set_zero_conjspinor
  interface operator (*)
     module procedure integer_spinor, spinor_integer, &
          real_spinor, double_spinor, &
          complex_spinor, dcomplex_spinor, &
          spinor_real, spinor_double, &
          spinor_complex, spinor_dcomplex
  end interface
  private :: integer_spinor, spinor_integer, real_spinor, &
       double_spinor, complex_spinor, dcomplex_spinor, &
       spinor_real, spinor_double, spinor_complex, spinor_dcomplex
  interface operator (*)
     module procedure integer_conjspinor, conjspinor_integer, &
          real_conjspinor, double_conjspinor, &
          complex_conjspinor, dcomplex_conjspinor, &
          conjspinor_real, conjspinor_double, &
          conjspinor_complex, conjspinor_dcomplex
  end interface
  private :: integer_conjspinor, conjspinor_integer, real_conjspinor, &
       double_conjspinor, complex_conjspinor, dcomplex_conjspinor, &
       conjspinor_real, conjspinor_double, conjspinor_complex, &
       conjspinor_dcomplex
  interface operator (+)
     module procedure plus_spinor, plus_conjspinor
  end interface
  private :: plus_spinor, plus_conjspinor
  interface operator (-)
     module procedure neg_spinor, neg_conjspinor
  end interface
  private :: neg_spinor, neg_conjspinor
  interface operator (+)
     module procedure add_spinor, add_conjspinor
  end interface
  private :: add_spinor, add_conjspinor
  interface operator (-)
     module procedure sub_spinor, sub_conjspinor
  end interface
  private :: sub_spinor, sub_conjspinor
  interface abs
     module procedure abs_spinor, abs_conjspinor
  end interface
  private :: abs_spinor, abs_conjspinor
  integer, parameter, public :: omega_spinors_2010_01_A = 0
contains
  pure function conjspinor_spinor (psibar, psi) result (psibarpsi)
    complex(kind=default) :: psibarpsi
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    psibarpsi = psibar%a(1)*psi%a(1) + psibar%a(2)*psi%a(2) &
              + psibar%a(3)*psi%a(3) + psibar%a(4)*psi%a(4)
  end function conjspinor_spinor
  elemental subroutine set_zero_spinor (x)
    type(spinor), intent(out) :: x
    x%a = 0
  end subroutine set_zero_spinor
  elemental subroutine set_zero_conjspinor (x)
    type(conjspinor), intent(out) :: x
    x%a = 0
  end subroutine set_zero_conjspinor
  pure function integer_spinor (x, y) result (xy)
    integer, intent(in) :: x
    type(spinor), intent(in) :: y
    type(spinor) :: xy
    xy%a = x * y%a
  end function integer_spinor
  pure function real_spinor (x, y) result (xy)
    real(kind=single), intent(in) :: x
    type(spinor), intent(in) :: y
    type(spinor) :: xy
    xy%a = x * y%a
  end function real_spinor
  pure function double_spinor (x, y) result (xy)
    real(kind=default), intent(in) :: x
    type(spinor), intent(in) :: y
    type(spinor) :: xy
    xy%a = x * y%a
  end function double_spinor
  pure function complex_spinor (x, y) result (xy)
    complex(kind=single), intent(in) :: x
    type(spinor), intent(in) :: y
    type(spinor) :: xy
    xy%a = x * y%a
  end function complex_spinor
  pure function dcomplex_spinor (x, y) result (xy)
    complex(kind=default), intent(in) :: x
    type(spinor), intent(in) :: y
    type(spinor) :: xy
    xy%a = x * y%a
  end function dcomplex_spinor
  pure function spinor_integer (y, x) result (xy)
    integer, intent(in) :: x
    type(spinor), intent(in) :: y
    type(spinor) :: xy
    xy%a = x * y%a
  end function spinor_integer
  pure function spinor_real (y, x) result (xy)
    real(kind=single), intent(in) :: x
    type(spinor), intent(in) :: y
    type(spinor) :: xy
    xy%a = x * y%a
  end function spinor_real
  pure function spinor_double (y, x) result (xy)
    real(kind=default), intent(in) :: x
    type(spinor), intent(in) :: y
    type(spinor) :: xy
    xy%a = x * y%a
  end function spinor_double
  pure function spinor_complex (y, x) result (xy)
    complex(kind=single), intent(in) :: x
    type(spinor), intent(in) :: y
    type(spinor) :: xy
    xy%a = x * y%a
  end function spinor_complex
  pure function spinor_dcomplex (y, x) result (xy)
    complex(kind=default), intent(in) :: x
    type(spinor), intent(in) :: y
    type(spinor) :: xy
    xy%a = x * y%a
  end function spinor_dcomplex
  pure function integer_conjspinor (x, y) result (xy)
    integer, intent(in) :: x
    type(conjspinor), intent(in) :: y
    type(conjspinor) :: xy
    xy%a = x * y%a
  end function integer_conjspinor
  pure function real_conjspinor (x, y) result (xy)
    real(kind=single), intent(in) :: x
    type(conjspinor), intent(in) :: y
    type(conjspinor) :: xy
    xy%a = x * y%a
  end function real_conjspinor
  pure function double_conjspinor (x, y) result (xy)
    real(kind=default), intent(in) :: x
    type(conjspinor), intent(in) :: y
    type(conjspinor) :: xy
    xy%a = x * y%a
  end function double_conjspinor
  pure function complex_conjspinor (x, y) result (xy)
    complex(kind=single), intent(in) :: x
    type(conjspinor), intent(in) :: y
    type(conjspinor) :: xy
    xy%a = x * y%a
  end function complex_conjspinor
  pure function dcomplex_conjspinor (x, y) result (xy)
    complex(kind=default), intent(in) :: x
    type(conjspinor), intent(in) :: y
    type(conjspinor) :: xy
    xy%a = x * y%a
  end function dcomplex_conjspinor
  pure function conjspinor_integer (y, x) result (xy)
    integer, intent(in) :: x
    type(conjspinor), intent(in) :: y
    type(conjspinor) :: xy
    xy%a = x * y%a
  end function conjspinor_integer
  pure function conjspinor_real (y, x) result (xy)
    real(kind=single), intent(in) :: x
    type(conjspinor), intent(in) :: y
    type(conjspinor) :: xy
    xy%a = x * y%a
  end function conjspinor_real
  pure function conjspinor_double (y, x) result (xy)
    real(kind=default), intent(in) :: x
    type(conjspinor), intent(in) :: y
    type(conjspinor) :: xy
    xy%a = x * y%a
  end function conjspinor_double
  pure function conjspinor_complex (y, x) result (xy)
    complex(kind=single), intent(in) :: x
    type(conjspinor), intent(in) :: y
    type(conjspinor) :: xy
    xy%a = x * y%a
  end function conjspinor_complex
  pure function conjspinor_dcomplex (y, x) result (xy)
    complex(kind=default), intent(in) :: x
    type(conjspinor), intent(in) :: y
    type(conjspinor) :: xy
    xy%a = x * y%a
  end function conjspinor_dcomplex
  pure function plus_spinor (x) result (plus_x)
    type(spinor), intent(in) :: x
    type(spinor) :: plus_x
    plus_x%a = x%a
  end function plus_spinor
  pure function neg_spinor (x) result (neg_x)
    type(spinor), intent(in) :: x
    type(spinor) :: neg_x
    neg_x%a = - x%a
  end function neg_spinor
  pure function plus_conjspinor (x) result (plus_x)
    type(conjspinor), intent(in) :: x
    type(conjspinor) :: plus_x
    plus_x%a = x%a
  end function plus_conjspinor
  pure function neg_conjspinor (x) result (neg_x)
    type(conjspinor), intent(in) :: x
    type(conjspinor) :: neg_x
    neg_x%a = - x%a
  end function neg_conjspinor
  pure function add_spinor (x, y) result (xy)
    type(spinor), intent(in) :: x, y
    type(spinor) :: xy
    xy%a = x%a + y%a
  end function add_spinor
  pure function sub_spinor (x, y) result (xy)
    type(spinor), intent(in) :: x, y
    type(spinor) :: xy
    xy%a = x%a - y%a
  end function sub_spinor
  pure function add_conjspinor (x, y) result (xy)
    type(conjspinor), intent(in) :: x, y
    type(conjspinor) :: xy
    xy%a = x%a + y%a
  end function add_conjspinor
  pure function sub_conjspinor (x, y) result (xy)
    type(conjspinor), intent(in) :: x, y
    type(conjspinor) :: xy
    xy%a = x%a - y%a
  end function sub_conjspinor
  pure function abs_spinor (psi) result (x)
    type(spinor), intent(in) :: psi
    real(kind=default) :: x
    x = sqrt (real (dot_product (psi%a, psi%a)))
  end function abs_spinor
  pure function abs_conjspinor (psibar) result (x)
    real(kind=default) :: x
    type(conjspinor), intent(in) :: psibar
    x = sqrt (real (dot_product (psibar%a, psibar%a)))
  end function abs_conjspinor
end module omega_spinors
