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
module omega_vectors
  use kinds
  use constants
  implicit none
  private
  public :: assignment (=), operator(==)
  public :: operator (*), operator (+), operator (-), operator (.wedge.)
  public :: abs, conjg, set_zero
  public :: random_momentum
  
  
  type, public :: momentum
     ! private (omegalib needs access, but DON'T TOUCH IT!)
     real(kind=default) :: t
     real(kind=default), dimension(3) :: x
  end type momentum
  type, public :: vector
     ! private (omegalib needs access, but DON'T TOUCH IT!)
     complex(kind=default) :: t
     complex(kind=default), dimension(3) :: x
  end type vector
  type, public :: tensor2odd
     ! private (omegalib needs access, but DON'T TOUCH IT!)
     complex(kind=default), dimension(3) :: e
     complex(kind=default), dimension(3) :: b
  end type tensor2odd
  interface assignment (=)
     module procedure momentum_of_array, vector_of_momentum, &
          vector_of_array, vector_of_double_array, &
          array_of_momentum, array_of_vector
  end interface
  private :: momentum_of_array, vector_of_momentum, vector_of_array, &
       vector_of_double_array, array_of_momentum, array_of_vector
  interface operator(==)
     module procedure momentum_eq
  end interface
  interface operator (*)
     module procedure momentum_momentum, vector_vector, &
          vector_momentum, momentum_vector, tensor2odd_tensor2odd
  end interface
  private :: momentum_momentum, vector_vector, vector_momentum, &
       momentum_vector, tensor2odd_tensor2odd
  interface operator (*)
     module procedure momentum_tensor2odd, tensor2odd_momentum, &
          vector_tensor2odd, tensor2odd_vector
  end interface
  private :: momentum_tensor2odd, tensor2odd_momentum, vector_tensor2odd, &
       tensor2odd_vector
  interface operator (.wedge.)
     module procedure momentum_wedge_momentum, &
          momentum_wedge_vector, vector_wedge_momentum, vector_wedge_vector
  end interface
  private :: momentum_wedge_momentum, momentum_wedge_vector, &
       vector_wedge_momentum, vector_wedge_vector
  interface set_zero
    module procedure set_zero_vector, set_zero_momentum, &
      set_zero_tensor2odd, set_zero_real, set_zero_complex
  end interface
  private :: set_zero_vector, set_zero_momentum, set_zero_tensor2odd
  interface operator (*)
     module procedure integer_momentum, real_momentum, double_momentum, &
          complex_momentum, dcomplex_momentum, &
          integer_vector, real_vector, double_vector, &
          complex_vector, dcomplex_vector, &
          integer_tensor2odd, real_tensor2odd, double_tensor2odd, &
          complex_tensor2odd, dcomplex_tensor2odd, &
          momentum_integer, momentum_real, momentum_double, &
          momentum_complex, momentum_dcomplex, &
          vector_integer, vector_real, vector_double, &
          vector_complex, vector_dcomplex, &
          tensor2odd_integer, tensor2odd_real, tensor2odd_double, &
          tensor2odd_complex, tensor2odd_dcomplex
  end interface
  private :: integer_momentum, real_momentum, double_momentum, &
       complex_momentum, dcomplex_momentum, integer_vector, real_vector, &
       double_vector, complex_vector, dcomplex_vector, &
       integer_tensor2odd, real_tensor2odd, double_tensor2odd, &
       complex_tensor2odd, dcomplex_tensor2odd, momentum_integer, &
       momentum_real, momentum_double, momentum_complex, &
       momentum_dcomplex, vector_integer, vector_real, vector_double, &
       vector_complex, vector_dcomplex, tensor2odd_integer, &
       tensor2odd_real, tensor2odd_double, tensor2odd_complex, &
       tensor2odd_dcomplex
  interface operator (+)
     module procedure plus_momentum, plus_vector, plus_tensor2odd
  end interface
  private :: plus_momentum, plus_vector, plus_tensor2odd
  interface operator (-)
     module procedure neg_momentum, neg_vector, neg_tensor2odd
  end interface
  private :: neg_momentum, neg_vector, neg_tensor2odd
  interface operator (+)
     module procedure add_momentum, add_vector, &
          add_vector_momentum, add_momentum_vector, add_tensor2odd
  end interface
  private :: add_momentum, add_vector, add_vector_momentum, &
       add_momentum_vector, add_tensor2odd
  interface operator (-)
     module procedure sub_momentum, sub_vector, &
          sub_vector_momentum, sub_momentum_vector, sub_tensor2odd
  end interface
  private :: sub_momentum, sub_vector, sub_vector_momentum, &
       sub_momentum_vector, sub_tensor2odd
  interface abs
     module procedure abs_momentum, abs_vector, abs_tensor2odd
  end interface
  private :: abs_momentum, abs_vector, abs_tensor2odd
  interface conjg
     module procedure conjg_momentum, conjg_vector, conjg_tensor2odd
  end interface
  private :: conjg_momentum, conjg_vector, conjg_tensor2odd
  interface pseudo_scalar
     module procedure pseudo_scalar_momentum, pseudo_scalar_vector, &
          pseudo_scalar_vec_mom
  end interface
  public :: pseudo_scalar
  private :: pseudo_scalar_momentum, pseudo_scalar_vector
  interface pseudo_vector
     module procedure pseudo_vector_momentum, pseudo_vector_vector, &
          pseudo_vector_vec_mom
  end interface
  public :: pseudo_vector
  private :: pseudo_vector_momentum, pseudo_vector_vector
  integer, parameter, public :: omega_vectors_2010_01_A = 0
contains
  pure subroutine momentum_of_array (m, p)
    type(momentum), intent(out) :: m
    real(kind=default), dimension(0:), intent(in) :: p
    m%t = p(0)
    m%x = p(1:3)
  end subroutine momentum_of_array
  pure subroutine array_of_momentum (p, v)
    real(kind=default), dimension(0:), intent(out) :: p
    type(momentum), intent(in) :: v
    p(0) = v%t
    p(1:3) = v%x
  end subroutine array_of_momentum
  pure subroutine vector_of_array (v, p)
    type(vector), intent(out) :: v
    complex(kind=default), dimension(0:), intent(in) :: p
    v%t = p(0)
    v%x = p(1:3)
  end subroutine vector_of_array
  pure subroutine vector_of_double_array (v, p)
    type(vector), intent(out) :: v
    real(kind=default), dimension(0:), intent(in) :: p
    v%t = p(0)
    v%x = p(1:3)
  end subroutine vector_of_double_array
  pure subroutine array_of_vector (p, v)
    complex(kind=default), dimension(0:), intent(out) :: p
    type(vector), intent(in) :: v
    p(0) = v%t
    p(1:3) = v%x
  end subroutine array_of_vector
  pure subroutine vector_of_momentum (v, p)
    type(vector), intent(out) :: v
    type(momentum), intent(in) :: p
    v%t = p%t
    v%x = p%x
  end subroutine vector_of_momentum
  elemental function momentum_eq (lhs, rhs) result (yorn)
    logical :: yorn
    type(momentum), intent(in) :: lhs
    type(momentum), intent(in) :: rhs
    yorn = all (abs(lhs%x - rhs%x) < eps0) .and. abs(lhs%t - rhs%t) < eps0
  end function momentum_eq
  pure function momentum_momentum (x, y) result (xy)
    type(momentum), intent(in) :: x
    type(momentum), intent(in) :: y
    real(kind=default) :: xy
    xy = x%t*y%t - x%x(1)*y%x(1) - x%x(2)*y%x(2) - x%x(3)*y%x(3)
  end function momentum_momentum
  pure function momentum_vector (x, y) result (xy)
    type(momentum), intent(in) :: x
    type(vector), intent(in) :: y
    complex(kind=default) :: xy
    xy = x%t*y%t - x%x(1)*y%x(1) - x%x(2)*y%x(2) - x%x(3)*y%x(3)
  end function momentum_vector
  pure function vector_momentum (x, y) result (xy)
    type(vector), intent(in) :: x
    type(momentum), intent(in) :: y
    complex(kind=default) :: xy
    xy = x%t*y%t - x%x(1)*y%x(1) - x%x(2)*y%x(2) - x%x(3)*y%x(3)
  end function vector_momentum
  pure function vector_vector (x, y) result (xy)
    type(vector), intent(in) :: x
    type(vector), intent(in) :: y
    complex(kind=default) :: xy
    xy = x%t*y%t - x%x(1)*y%x(1) - x%x(2)*y%x(2) - x%x(3)*y%x(3)
  end function vector_vector
  pure function tensor2odd_tensor2odd (x, y) result (xy)
    type(tensor2odd), intent(in) :: x
    type(tensor2odd), intent(in) :: y
    complex(kind=default) :: xy
    xy = x%b(1)*y%b(1) + x%b(2)*y%b(2) + x%b(3)*y%b(3) &
       - x%e(1)*y%e(1) - x%e(2)*y%e(2) - x%e(3)*y%e(3)
  end function tensor2odd_tensor2odd
  pure function vector_tensor2odd (x, t2) result (xt2)
    type(vector), intent(in) :: x
    type(tensor2odd), intent(in) :: t2
    type(vector) :: xt2
    xt2%t = x%x(1)*t2%e(1) + x%x(2)*t2%e(2) + x%x(3)*t2%e(3)
    xt2%x(1) = x%t*t2%e(1) + x%x(2)*t2%b(3) - x%x(3)*t2%b(2)
    xt2%x(2) = x%t*t2%e(2) + x%x(3)*t2%b(1) - x%x(1)*t2%b(3)
    xt2%x(3) = x%t*t2%e(3) + x%x(1)*t2%b(2) - x%x(2)*t2%b(1)
  end function vector_tensor2odd
  pure function momentum_tensor2odd (x, t2) result (xt2)
    type(momentum), intent(in) :: x
    type(tensor2odd), intent(in) :: t2
    type(vector) :: xt2
    xt2%t = x%x(1)*t2%e(1) + x%x(2)*t2%e(2) + x%x(3)*t2%e(3)
    xt2%x(1) = x%t*t2%e(1) + x%x(2)*t2%b(3) - x%x(3)*t2%b(2)
    xt2%x(2) = x%t*t2%e(2) + x%x(3)*t2%b(1) - x%x(1)*t2%b(3)
    xt2%x(3) = x%t*t2%e(3) + x%x(1)*t2%b(2) - x%x(2)*t2%b(1)
  end function momentum_tensor2odd
  pure function tensor2odd_vector (t2, x) result (t2x)
    type(tensor2odd), intent(in) :: t2
    type(vector), intent(in) :: x
    type(vector) :: t2x
    t2x%t = - t2%e(1)*x%x(1) - t2%e(2)*x%x(2) - t2%e(3)*x%x(3)
    t2x%x(1) = - t2%e(1)*x%t + t2%b(2)*x%x(3) - t2%b(3)*x%x(2)
    t2x%x(2) = - t2%e(2)*x%t + t2%b(3)*x%x(1) - t2%b(1)*x%x(3)
    t2x%x(3) = - t2%e(3)*x%t + t2%b(1)*x%x(2) - t2%b(2)*x%x(1)
  end function tensor2odd_vector
  pure function tensor2odd_momentum (t2, x) result (t2x)
    type(tensor2odd), intent(in) :: t2
    type(momentum), intent(in) :: x
    type(vector) :: t2x
    t2x%t = - t2%e(1)*x%x(1) - t2%e(2)*x%x(2) - t2%e(3)*x%x(3)
    t2x%x(1) = - t2%e(1)*x%t + t2%b(2)*x%x(3) - t2%b(3)*x%x(2)
    t2x%x(2) = - t2%e(2)*x%t + t2%b(3)*x%x(1) - t2%b(1)*x%x(3)
    t2x%x(3) = - t2%e(3)*x%t + t2%b(1)*x%x(2) - t2%b(2)*x%x(1)
  end function tensor2odd_momentum
  pure function momentum_wedge_momentum (x, y) result (t2)
    type(momentum), intent(in) :: x
    type(momentum), intent(in) :: y
    type(tensor2odd) :: t2
    t2%e = x%t * y%x - x%x * y%t
    t2%b(1) = x%x(2) * y%x(3) - x%x(3) * y%x(2)
    t2%b(2) = x%x(3) * y%x(1) - x%x(1) * y%x(3)
    t2%b(3) = x%x(1) * y%x(2) - x%x(2) * y%x(1)
  end function momentum_wedge_momentum
  pure function momentum_wedge_vector (x, y) result (t2)
    type(momentum), intent(in) :: x
    type(vector), intent(in) :: y
    type(tensor2odd) :: t2
    t2%e = x%t * y%x - x%x * y%t
    t2%b(1) = x%x(2) * y%x(3) - x%x(3) * y%x(2)
    t2%b(2) = x%x(3) * y%x(1) - x%x(1) * y%x(3)
    t2%b(3) = x%x(1) * y%x(2) - x%x(2) * y%x(1)
  end function momentum_wedge_vector
  pure function vector_wedge_momentum (x, y) result (t2)
    type(vector), intent(in) :: x
    type(momentum), intent(in) :: y
    type(tensor2odd) :: t2
    t2%e = x%t * y%x - x%x * y%t
    t2%b(1) = x%x(2) * y%x(3) - x%x(3) * y%x(2)
    t2%b(2) = x%x(3) * y%x(1) - x%x(1) * y%x(3)
    t2%b(3) = x%x(1) * y%x(2) - x%x(2) * y%x(1)
  end function vector_wedge_momentum
  pure function vector_wedge_vector (x, y) result (t2)
    type(vector), intent(in) :: x
    type(vector), intent(in) :: y
    type(tensor2odd) :: t2
    t2%e = x%t * y%x - x%x * y%t
    t2%b(1) = x%x(2) * y%x(3) - x%x(3) * y%x(2)
    t2%b(2) = x%x(3) * y%x(1) - x%x(1) * y%x(3)
    t2%b(3) = x%x(1) * y%x(2) - x%x(2) * y%x(1)
  end function vector_wedge_vector
  elemental subroutine set_zero_vector (x)
    type(vector), intent(out) :: x
    x%t = 0
    x%x = 0
  end subroutine set_zero_vector
  elemental subroutine set_zero_momentum (x)
    type(momentum), intent(out) :: x
    x%t = 0
    x%x = 0
  end subroutine set_zero_momentum
  elemental subroutine set_zero_tensor2odd (x)
    type(tensor2odd), intent(out) :: x
    x%e = 0
    x%b = 0
  end subroutine set_zero_tensor2odd
  elemental subroutine set_zero_real (x)
    real(kind=default), intent(out) :: x
    x = 0
  end subroutine set_zero_real
  elemental subroutine set_zero_complex (x)
    complex(kind=default), intent(out) :: x
    x = 0
  end subroutine set_zero_complex
  pure function integer_momentum (x, y) result (xy)
    integer, intent(in) :: x
    type(momentum), intent(in) :: y
    type(momentum) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function integer_momentum
  pure function real_momentum (x, y) result (xy)
    real(kind=single), intent(in) :: x
    type(momentum), intent(in) :: y
    type(momentum) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function real_momentum
  pure function double_momentum (x, y) result (xy)
    real(kind=default), intent(in) :: x
    type(momentum), intent(in) :: y
    type(momentum) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function double_momentum
  pure function complex_momentum (x, y) result (xy)
    complex(kind=single), intent(in) :: x
    type(momentum), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function complex_momentum
  pure function dcomplex_momentum (x, y) result (xy)
    complex(kind=default), intent(in) :: x
    type(momentum), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function dcomplex_momentum
  pure function integer_vector (x, y) result (xy)
    integer, intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function integer_vector
  pure function real_vector (x, y) result (xy)
    real(kind=single), intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function real_vector
  pure function double_vector (x, y) result (xy)
    real(kind=default), intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function double_vector
  pure function complex_vector (x, y) result (xy)
    complex(kind=single), intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function complex_vector
  pure function dcomplex_vector (x, y) result (xy)
    complex(kind=default), intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function dcomplex_vector
  pure function integer_tensor2odd (x, t2) result (xt2)
    integer, intent(in) :: x
    type(tensor2odd), intent(in) :: t2
    type(tensor2odd) :: xt2
    xt2%e = x * t2%e
    xt2%b = x * t2%b
  end function integer_tensor2odd
  pure function real_tensor2odd (x, t2) result (xt2)
    real(kind=single), intent(in) :: x
    type(tensor2odd), intent(in) :: t2
    type(tensor2odd) :: xt2
    xt2%e = x * t2%e
    xt2%b = x * t2%b
  end function real_tensor2odd
  pure function double_tensor2odd (x, t2) result (xt2)
    real(kind=default), intent(in) :: x
    type(tensor2odd), intent(in) :: t2
    type(tensor2odd) :: xt2
    xt2%e = x * t2%e
    xt2%b = x * t2%b
  end function double_tensor2odd
  pure function complex_tensor2odd (x, t2) result (xt2)
    complex(kind=single), intent(in) :: x
    type(tensor2odd), intent(in) :: t2
    type(tensor2odd) :: xt2
    xt2%e = x * t2%e
    xt2%b = x * t2%b
  end function complex_tensor2odd
  pure function dcomplex_tensor2odd (x, t2) result (xt2)
    complex(kind=default), intent(in) :: x
    type(tensor2odd), intent(in) :: t2
    type(tensor2odd) :: xt2
    xt2%e = x * t2%e
    xt2%b = x * t2%b
  end function dcomplex_tensor2odd
  pure function momentum_integer (y, x) result (xy)
    integer, intent(in) :: x
    type(momentum), intent(in) :: y
    type(momentum) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function momentum_integer
  pure function momentum_real (y, x) result (xy)
    real(kind=single), intent(in) :: x
    type(momentum), intent(in) :: y
    type(momentum) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function momentum_real
  pure function momentum_double (y, x) result (xy)
    real(kind=default), intent(in) :: x
    type(momentum), intent(in) :: y
    type(momentum) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function momentum_double
  pure function momentum_complex (y, x) result (xy)
    complex(kind=single), intent(in) :: x
    type(momentum), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function momentum_complex
  pure function momentum_dcomplex (y, x) result (xy)
    complex(kind=default), intent(in) :: x
    type(momentum), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function momentum_dcomplex
  pure function vector_integer (y, x) result (xy)
    integer, intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function vector_integer
  pure function vector_real (y, x) result (xy)
    real(kind=single), intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function vector_real
  pure function vector_double (y, x) result (xy)
    real(kind=default), intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function vector_double
  pure function vector_complex (y, x) result (xy)
    complex(kind=single), intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function vector_complex
  pure function vector_dcomplex (y, x) result (xy)
    complex(kind=default), intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x * y%t
    xy%x = x * y%x
  end function vector_dcomplex
  pure function tensor2odd_integer (t2, x) result (t2x)
    type(tensor2odd), intent(in) :: t2
    integer, intent(in) :: x
    type(tensor2odd) :: t2x
    t2x%e = x * t2%e
    t2x%b = x * t2%b
  end function tensor2odd_integer
  pure function tensor2odd_real (t2, x) result (t2x)
    type(tensor2odd), intent(in) :: t2
    real(kind=single), intent(in) :: x
    type(tensor2odd) :: t2x
    t2x%e = x * t2%e
    t2x%b = x * t2%b
  end function tensor2odd_real
  pure function tensor2odd_double (t2, x) result (t2x)
    type(tensor2odd), intent(in) :: t2
    real(kind=default), intent(in) :: x
    type(tensor2odd) :: t2x
    t2x%e = x * t2%e
    t2x%b = x * t2%b
  end function tensor2odd_double
  pure function tensor2odd_complex (t2, x) result (t2x)
    type(tensor2odd), intent(in) :: t2
    complex(kind=single), intent(in) :: x
    type(tensor2odd) :: t2x
    t2x%e = x * t2%e
    t2x%b = x * t2%b
  end function tensor2odd_complex
  pure function tensor2odd_dcomplex (t2, x) result (t2x)
    type(tensor2odd), intent(in) :: t2
    complex(kind=default), intent(in) :: x
    type(tensor2odd) :: t2x
    t2x%e = x * t2%e
    t2x%b = x * t2%b
  end function tensor2odd_dcomplex
  pure function plus_momentum (x) result (plus_x)
    type(momentum), intent(in) :: x
    type(momentum) :: plus_x
    plus_x = x
  end function plus_momentum
  pure function neg_momentum (x) result (neg_x)
    type(momentum), intent(in) :: x
    type(momentum) :: neg_x
    neg_x%t = - x%t
    neg_x%x = - x%x
  end function neg_momentum
  pure function plus_vector (x) result (plus_x)
    type(vector), intent(in) :: x
    type(vector) :: plus_x
    plus_x = x
  end function plus_vector
  pure function neg_vector (x) result (neg_x)
    type(vector), intent(in) :: x
    type(vector) :: neg_x
    neg_x%t = - x%t
    neg_x%x = - x%x
  end function neg_vector
  pure function plus_tensor2odd (x) result (plus_x)
    type(tensor2odd), intent(in) :: x
    type(tensor2odd) :: plus_x
    plus_x = x
  end function plus_tensor2odd
  pure function neg_tensor2odd (x) result (neg_x)
    type(tensor2odd), intent(in) :: x
    type(tensor2odd) :: neg_x
    neg_x%e = - x%e
    neg_x%b = - x%b
  end function neg_tensor2odd
  pure function add_momentum (x, y) result (xy)
    type(momentum), intent(in) :: x, y
    type(momentum) :: xy
    xy%t = x%t + y%t
    xy%x = x%x + y%x
  end function add_momentum
  pure function add_vector (x, y) result (xy)
    type(vector), intent(in) :: x, y
    type(vector) :: xy
    xy%t = x%t + y%t
    xy%x = x%x + y%x
  end function add_vector
  pure function add_momentum_vector (x, y) result (xy)
    type(momentum), intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x%t + y%t
    xy%x = x%x + y%x
  end function add_momentum_vector
  pure function add_vector_momentum (x, y) result (xy)
    type(vector), intent(in) :: x
    type(momentum), intent(in) :: y
    type(vector) :: xy
    xy%t = x%t + y%t
    xy%x = x%x + y%x
  end function add_vector_momentum
  pure function add_tensor2odd (x, y) result (xy)
    type(tensor2odd), intent(in) :: x, y
    type(tensor2odd) :: xy
    xy%e = x%e + y%e
    xy%b = x%b + y%b
  end function add_tensor2odd
  pure function sub_momentum (x, y) result (xy)
    type(momentum), intent(in) :: x, y
    type(momentum) :: xy
    xy%t = x%t - y%t
    xy%x = x%x - y%x
  end function sub_momentum
  pure function sub_vector (x, y) result (xy)
    type(vector), intent(in) :: x, y
    type(vector) :: xy
    xy%t = x%t - y%t
    xy%x = x%x - y%x
  end function sub_vector
  pure function sub_momentum_vector (x, y) result (xy)
    type(momentum), intent(in) :: x
    type(vector), intent(in) :: y
    type(vector) :: xy
    xy%t = x%t - y%t
    xy%x = x%x - y%x
  end function sub_momentum_vector
  pure function sub_vector_momentum (x, y) result (xy)
    type(vector), intent(in) :: x
    type(momentum), intent(in) :: y
    type(vector) :: xy
    xy%t = x%t - y%t
    xy%x = x%x - y%x
  end function sub_vector_momentum
  pure function sub_tensor2odd (x, y) result (xy)
    type(tensor2odd), intent(in) :: x, y
    type(tensor2odd) :: xy
    xy%e = x%e - y%e
    xy%b = x%b - y%b
  end function sub_tensor2odd
  pure function abs_momentum (x) result (absx)
    type(momentum), intent(in) :: x
    real(kind=default) :: absx
    absx = sqrt (real (x%t*x%t + dot_product (x%x, x%x)))
  end function abs_momentum
  pure function abs_vector (x) result (absx)
    type(vector), intent(in) :: x
    real(kind=default) :: absx
    absx = sqrt (real (conjg(x%t)*x%t + dot_product (x%x, x%x)))
  end function abs_vector
  pure function abs_tensor2odd (x) result (absx)
    type(tensor2odd), intent(in) :: x
    real(kind=default) :: absx
    absx = sqrt (real (dot_product (x%e, x%e) + dot_product (x%b, x%b)))
  end function abs_tensor2odd
  pure function conjg_momentum (x) result (conjg_x)
    type(momentum), intent(in) :: x
    type(momentum) :: conjg_x
    conjg_x = x
  end function conjg_momentum
  pure function conjg_vector (x) result (conjg_x)
    type(vector), intent(in) :: x
    type(vector) :: conjg_x
    conjg_x%t = conjg (x%t)
    conjg_x%x = conjg (x%x)
  end function conjg_vector
  pure function conjg_tensor2odd (t2) result (conjg_t2)
    type(tensor2odd), intent(in) :: t2
    type(tensor2odd) :: conjg_t2
    conjg_t2%e = conjg (t2%e)
    conjg_t2%b = conjg (t2%b)
  end function conjg_tensor2odd
  pure function pseudo_scalar_momentum (p1, p2, p3, p4) result (eps1234)
    type(momentum), intent(in) :: p1, p2, p3, p4
    real(kind=default) :: eps1234
    eps1234 = &
         p1%t    * p2%x(1) * (p3%x(2) * p4%x(3) - p3%x(3) * p4%x(2)) &
       + p1%t    * p2%x(2) * (p3%x(3) * p4%x(1) - p3%x(1) * p4%x(3)) &
       + p1%t    * p2%x(3) * (p3%x(1) * p4%x(2) - p3%x(2) * p4%x(1)) &
       - p1%x(1) * p2%x(2) * (p3%x(3) * p4%t    - p3%t    * p4%x(3)) &
       - p1%x(1) * p2%x(3) * (p3%t    * p4%x(2) - p3%x(2) * p4%t   ) &
       - p1%x(1) * p2%t    * (p3%x(2) * p4%x(3) - p3%x(3) * p4%x(2)) &
       + p1%x(2) * p2%x(3) * (p3%t    * p4%x(1) - p3%x(1) * p4%t   ) &
       + p1%x(2) * p2%t    * (p3%x(1) * p4%x(3) - p3%x(3) * p4%x(1)) &
       + p1%x(2) * p2%x(1) * (p3%x(3) * p4%t    - p3%t    * p4%x(3)) &
       - p1%x(3) * p2%t    * (p3%x(1) * p4%x(2) - p3%x(2) * p4%x(1)) &
       - p1%x(3) * p2%x(1) * (p3%x(2) * p4%t    - p3%t    * p4%x(2)) &
       - p1%x(3) * p2%x(2) * (p3%t    * p4%x(1) - p3%x(1) * p4%t   )
  end function pseudo_scalar_momentum
  pure function pseudo_scalar_vector (p1, p2, p3, p4) result (eps1234)
    type(vector), intent(in) :: p1, p2, p3, p4
    complex(kind=default) :: eps1234
    eps1234 = &
         p1%t    * p2%x(1) * (p3%x(2) * p4%x(3) - p3%x(3) * p4%x(2)) &
       + p1%t    * p2%x(2) * (p3%x(3) * p4%x(1) - p3%x(1) * p4%x(3)) &
       + p1%t    * p2%x(3) * (p3%x(1) * p4%x(2) - p3%x(2) * p4%x(1)) &
       - p1%x(1) * p2%x(2) * (p3%x(3) * p4%t    - p3%t    * p4%x(3)) &
       - p1%x(1) * p2%x(3) * (p3%t    * p4%x(2) - p3%x(2) * p4%t   ) &
       - p1%x(1) * p2%t    * (p3%x(2) * p4%x(3) - p3%x(3) * p4%x(2)) &
       + p1%x(2) * p2%x(3) * (p3%t    * p4%x(1) - p3%x(1) * p4%t   ) &
       + p1%x(2) * p2%t    * (p3%x(1) * p4%x(3) - p3%x(3) * p4%x(1)) &
       + p1%x(2) * p2%x(1) * (p3%x(3) * p4%t    - p3%t    * p4%x(3)) &
       - p1%x(3) * p2%t    * (p3%x(1) * p4%x(2) - p3%x(2) * p4%x(1)) &
       - p1%x(3) * p2%x(1) * (p3%x(2) * p4%t    - p3%t    * p4%x(2)) &
       - p1%x(3) * p2%x(2) * (p3%t    * p4%x(1) - p3%x(1) * p4%t   )
  end function pseudo_scalar_vector
  pure function pseudo_scalar_vec_mom (p1, v1, p2, v2) result (eps1234)
    type(momentum), intent(in)   :: p1, p2
    type(vector), intent(in) :: v1, v2
    complex(kind=default) :: eps1234
    eps1234 = &
         p1%t    * v1%x(1) * (p2%x(2) * v2%x(3) - p2%x(3) * v2%x(2)) &
       + p1%t    * v1%x(2) * (p2%x(3) * v2%x(1) - p2%x(1) * v2%x(3)) &
       + p1%t    * v1%x(3) * (p2%x(1) * v2%x(2) - p2%x(2) * v2%x(1)) &
       - p1%x(1) * v1%x(2) * (p2%x(3) * v2%t    - p2%t    * v2%x(3)) &
       - p1%x(1) * v1%x(3) * (p2%t    * v2%x(2) - p2%x(2) * v2%t   ) &
       - p1%x(1) * v1%t    * (p2%x(2) * v2%x(3) - p2%x(3) * v2%x(2)) &
       + p1%x(2) * v1%x(3) * (p2%t    * v2%x(1) - p2%x(1) * v2%t   ) &
       + p1%x(2) * v1%t    * (p2%x(1) * v2%x(3) - p2%x(3) * v2%x(1)) &
       + p1%x(2) * v1%x(1) * (p2%x(3) * v2%t    - p2%t    * v2%x(3)) &
       - p1%x(3) * v1%t    * (p2%x(1) * v2%x(2) - p2%x(2) * v2%x(1)) &
       - p1%x(3) * v1%x(1) * (p2%x(2) * v2%t    - p2%t    * v2%x(2)) &
       - p1%x(3) * v1%x(2) * (p2%t    * v2%x(1) - p2%x(1) * v2%t   )
  end function pseudo_scalar_vec_mom
  pure function pseudo_vector_momentum (p1, p2, p3) result (eps123)
    type(momentum), intent(in) :: p1, p2, p3
    type(momentum) :: eps123
    eps123%t = &
      + p1%x(1) * (p2%x(2) * p3%x(3) - p2%x(3) * p3%x(2)) &
      + p1%x(2) * (p2%x(3) * p3%x(1) - p2%x(1) * p3%x(3)) &
      + p1%x(3) * (p2%x(1) * p3%x(2) - p2%x(2) * p3%x(1))
    eps123%x(1) = &
      + p1%x(2) * (p2%x(3) * p3%t    - p2%t    * p3%x(3)) &
      + p1%x(3) * (p2%t    * p3%x(2) - p2%x(2) * p3%t   ) &
      + p1%t    * (p2%x(2) * p3%x(3) - p2%x(3) * p3%x(2))
    eps123%x(2) = &
      - p1%x(3) * (p2%t    * p3%x(1) - p2%x(1) * p3%t   ) &
      - p1%t    * (p2%x(1) * p3%x(3) - p2%x(3) * p3%x(1)) &
      - p1%x(1) * (p2%x(3) * p3%t    - p2%t    * p3%x(3))
    eps123%x(3) =  &
      + p1%t    * (p2%x(1) * p3%x(2) - p2%x(2) * p3%x(1)) &
      + p1%x(1) * (p2%x(2) * p3%t    - p2%t    * p3%x(2)) &
      + p1%x(2) * (p2%t    * p3%x(1) - p2%x(1) * p3%t   )
  end function pseudo_vector_momentum
  pure function pseudo_vector_vector (p1, p2, p3) result (eps123)
    type(vector), intent(in) :: p1, p2, p3
    type(vector) :: eps123
    eps123%t = &
      + p1%x(1) * (p2%x(2) * p3%x(3) - p2%x(3) * p3%x(2)) &
      + p1%x(2) * (p2%x(3) * p3%x(1) - p2%x(1) * p3%x(3)) &
      + p1%x(3) * (p2%x(1) * p3%x(2) - p2%x(2) * p3%x(1))
    eps123%x(1) = &
      + p1%x(2) * (p2%x(3) * p3%t    - p2%t    * p3%x(3)) &
      + p1%x(3) * (p2%t    * p3%x(2) - p2%x(2) * p3%t   ) &
      + p1%t    * (p2%x(2) * p3%x(3) - p2%x(3) * p3%x(2))
    eps123%x(2) = &
      - p1%x(3) * (p2%t    * p3%x(1) - p2%x(1) * p3%t   ) &
      - p1%t    * (p2%x(1) * p3%x(3) - p2%x(3) * p3%x(1)) &
      - p1%x(1) * (p2%x(3) * p3%t    - p2%t    * p3%x(3))
    eps123%x(3) =  &
      + p1%t    * (p2%x(1) * p3%x(2) - p2%x(2) * p3%x(1)) &
      + p1%x(1) * (p2%x(2) * p3%t    - p2%t    * p3%x(2)) &
      + p1%x(2) * (p2%t    * p3%x(1) - p2%x(1) * p3%t   )
  end function pseudo_vector_vector
  pure function pseudo_vector_vec_mom (p1, p2, v) result (eps123)
    type(momentum), intent(in) :: p1, p2
    type(vector), intent(in)   :: v
    type(vector) :: eps123
    eps123%t = &
      + p1%x(1) * (p2%x(2) * v%x(3) - p2%x(3) * v%x(2)) &
      + p1%x(2) * (p2%x(3) * v%x(1) - p2%x(1) * v%x(3)) &
      + p1%x(3) * (p2%x(1) * v%x(2) - p2%x(2) * v%x(1))
    eps123%x(1) = &
      + p1%x(2) * (p2%x(3) * v%t    - p2%t    * v%x(3)) &
      + p1%x(3) * (p2%t    * v%x(2) - p2%x(2) * v%t   ) &
      + p1%t    * (p2%x(2) * v%x(3) - p2%x(3) * v%x(2))
    eps123%x(2) = &
      - p1%x(3) * (p2%t    * v%x(1) - p2%x(1) * v%t   ) &
      - p1%t    * (p2%x(1) * v%x(3) - p2%x(3) * v%x(1)) &
      - p1%x(1) * (p2%x(3) * v%t    - p2%t    * v%x(3))
    eps123%x(3) =  &
      + p1%t    * (p2%x(1) * v%x(2) - p2%x(2) * v%x(1)) &
      + p1%x(1) * (p2%x(2) * v%t    - p2%t    * v%x(2)) &
      + p1%x(2) * (p2%t    * v%x(1) - p2%x(1) * v%t   )
  end function pseudo_vector_vec_mom
  subroutine random_momentum (p, pabs, m)
    type(momentum), intent(out) :: p
    real(kind=default), intent(in) :: pabs, m
    real(kind=default), dimension(2) :: r
    real(kind=default) :: phi, cos_th
    call random_number (r)
    phi = 2*PI * r(1)
    cos_th = 2 * r(2) - 1
    p%t = sqrt (pabs**2 + m**2)
    p%x = pabs * (/ cos_th * cos(phi), cos_th * sin(phi), sqrt (1 - cos_th**2) /)
  end subroutine random_momentum
end module omega_vectors
