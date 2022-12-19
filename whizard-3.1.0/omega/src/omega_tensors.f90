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
module omega_tensors
  use kinds
  use constants
  use omega_vectors
  implicit none
  private
  public :: operator (*), operator (+), operator (-), &
       operator (.tprod.)
  public :: abs, conjg, set_zero
  
  
  type, public :: tensor
     ! private (omegalib needs access, but DON'T TOUCH IT!)
     complex(kind=default), dimension(0:3,0:3) :: t
  end type tensor
  interface set_zero
    module procedure set_zero_tensor
  end interface
  private :: set_zero_tensor
  interface operator (*)
     module procedure integer_tensor, real_tensor, double_tensor, &
          complex_tensor, dcomplex_tensor
  end interface
  private :: integer_tensor, real_tensor, double_tensor
  private :: complex_tensor, dcomplex_tensor
  interface operator (+)
     module procedure plus_tensor
  end interface
  private :: plus_tensor
  interface operator (-)
    module procedure neg_tensor
  end interface
  private :: neg_tensor
  interface operator (+)
     module procedure add_tensor
  end interface
  private :: add_tensor
  interface operator (-)
     module procedure sub_tensor
  end interface
  private :: sub_tensor
  interface operator (.tprod.)
     module procedure out_prod_vv, out_prod_vm, &
          out_prod_mv, out_prod_mm
  end interface
  private :: out_prod_vv, out_prod_vm, &
       out_prod_mv, out_prod_mm
  interface abs
     module procedure abs_tensor
  end interface
  private :: abs_tensor
  interface conjg
     module procedure conjg_tensor
  end interface
  private :: conjg_tensor
  interface operator (*)
     module procedure tensor_tensor, vector_tensor, tensor_vector, &
          momentum_tensor, tensor_momentum
  end interface
  private :: tensor_tensor, vector_tensor, tensor_vector, &
       momentum_tensor, tensor_momentum
  integer, parameter, public :: omega_tensors_2010_01_A = 0
contains
  elemental subroutine set_zero_tensor (x)
    type(tensor), intent(out) :: x
    x%t = 0
  end subroutine set_zero_tensor
  pure function integer_tensor (x, y) result (xy)
    integer, intent(in) :: x
    type(tensor), intent(in) :: y
    type(tensor) :: xy
    xy%t = x * y%t
  end function integer_tensor
  pure function real_tensor (x, y) result (xy)
    real(kind=single), intent(in) :: x
    type(tensor), intent(in) :: y
    type(tensor) :: xy
    xy%t = x * y%t
  end function real_tensor
  pure function double_tensor (x, y) result (xy)
    real(kind=default), intent(in) :: x
    type(tensor), intent(in) :: y
    type(tensor) :: xy
    xy%t = x * y%t
  end function double_tensor
  pure function complex_tensor (x, y) result (xy)
    complex(kind=single), intent(in) :: x
    type(tensor), intent(in) :: y
    type(tensor) :: xy
    xy%t = x * y%t
  end function complex_tensor
  pure function dcomplex_tensor (x, y) result (xy)
    complex(kind=default), intent(in) :: x
    type(tensor), intent(in) :: y
    type(tensor) :: xy
    xy%t = x * y%t
  end function dcomplex_tensor
  pure function plus_tensor (t1) result (t2)
    type(tensor), intent(in) :: t1
    type(tensor) :: t2
    t2 = t1
  end function plus_tensor
  pure function neg_tensor (t1) result (t2)
    type(tensor), intent(in) :: t1
    type(tensor) :: t2
    t2%t = - t1%t
  end function neg_tensor
  pure function add_tensor (x, y) result (xy)
    type(tensor), intent(in) :: x, y
    type(tensor) :: xy
    xy%t = x%t + y%t
  end function add_tensor
  pure function sub_tensor (x, y) result (xy)
    type(tensor), intent(in) :: x, y
    type(tensor) :: xy
    xy%t = x%t - y%t
  end function sub_tensor
  pure function out_prod_vv (v, w) result (t)
    type(tensor) :: t
    type(vector), intent(in) :: v, w
    integer :: i, j
    t%t(0,0) = v%t * w%t
    t%t(0,1:3) = v%t * w%x
    t%t(1:3,0) = v%x * w%t
    do i = 1, 3
       do j = 1, 3
          t%t(i,j) = v%x(i) * w%x(j)
       end do
    end do
  end function out_prod_vv
  pure function out_prod_vm (v, m) result (t)
    type(tensor) :: t
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: m
    integer :: i, j
    t%t(0,0) = v%t * m%t
    t%t(0,1:3) = v%t * m%x
    t%t(1:3,0) = v%x * m%t
    do i = 1, 3
       do j = 1, 3
          t%t(i,j) = v%x(i) * m%x(j)
       end do
    end do
  end function out_prod_vm
  pure function out_prod_mv (m, v) result (t)
    type(tensor) :: t
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: m
    integer :: i, j
    t%t(0,0) = m%t * v%t
    t%t(0,1:3) = m%t * v%x
    t%t(1:3,0) = m%x * v%t
    do i = 1, 3
       do j = 1, 3
          t%t(i,j) = m%x(i) * v%x(j)
       end do
    end do
  end function out_prod_mv
  pure function out_prod_mm (m, n) result (t)
    type(tensor) :: t
    type(momentum), intent(in) :: m, n
    integer :: i, j
    t%t(0,0) = m%t * n%t
    t%t(0,1:3) = m%t * n%x
    t%t(1:3,0) = m%x * n%t
    do i = 1, 3
       do j = 1, 3
          t%t(i,j) = m%x(i) * n%x(j)
       end do
    end do
  end function out_prod_mm
  pure function abs_tensor (t) result (abs_t)
    type(tensor), intent(in) :: t
    real(kind=default) :: abs_t
    abs_t = sqrt (sum ((abs (t%t))**2))
  end function abs_tensor
  pure function conjg_tensor (t) result (conjg_t)
    type(tensor), intent(in) :: t
    type(tensor) :: conjg_t
    conjg_t%t = conjg (t%t)
  end function conjg_tensor
  pure function tensor_tensor (t1, t2) result (t1t2)
    type(tensor), intent(in) :: t1
    type(tensor), intent(in) :: t2
    complex(kind=default) :: t1t2
    integer :: i1, i2
    t1t2 = t1%t(0,0)*t2%t(0,0) &
         - dot_product (conjg (t1%t(0,1:)), t2%t(0,1:)) &
         - dot_product (conjg (t1%t(1:,0)), t2%t(1:,0))
    do i1 = 1, 3
       do i2 = 1, 3
          t1t2 = t1t2 + t1%t(i1,i2)*t2%t(i1,i2)
       end do
    end do
  end function tensor_tensor
  pure function tensor_vector (t, v) result (tv)
    type(tensor), intent(in) :: t
    type(vector), intent(in) :: v
    type(vector) :: tv
    tv%t =    t%t(0,0) * v%t - dot_product (conjg (t%t(0,1:)), v%x)
    tv%x(1) = t%t(0,1) * v%t - dot_product (conjg (t%t(1,1:)), v%x)
    tv%x(2) = t%t(0,2) * v%t - dot_product (conjg (t%t(2,1:)), v%x)
    tv%x(3) = t%t(0,3) * v%t - dot_product (conjg (t%t(3,1:)), v%x)
  end function tensor_vector
  pure function vector_tensor (v, t) result (vt)
    type(vector), intent(in) :: v
    type(tensor), intent(in) :: t
    type(vector) :: vt
    vt%t =    v%t * t%t(0,0) - dot_product (conjg (v%x), t%t(1:,0))
    vt%x(1) = v%t * t%t(0,1) - dot_product (conjg (v%x), t%t(1:,1))
    vt%x(2) = v%t * t%t(0,2) - dot_product (conjg (v%x), t%t(1:,2))
    vt%x(3) = v%t * t%t(0,3) - dot_product (conjg (v%x), t%t(1:,3))
  end function vector_tensor
  pure function tensor_momentum (t, p) result (tp)
    type(tensor), intent(in) :: t
    type(momentum), intent(in) :: p
    type(vector) :: tp
    tp%t =    t%t(0,0) * p%t - dot_product (conjg (t%t(0,1:)), p%x)
    tp%x(1) = t%t(0,1) * p%t - dot_product (conjg (t%t(1,1:)), p%x)
    tp%x(2) = t%t(0,2) * p%t - dot_product (conjg (t%t(2,1:)), p%x)
    tp%x(3) = t%t(0,3) * p%t - dot_product (conjg (t%t(3,1:)), p%x)
  end function tensor_momentum
  pure function momentum_tensor (p, t) result (pt)
    type(momentum), intent(in) :: p
    type(tensor), intent(in) :: t
    type(vector) :: pt
    pt%t =    p%t * t%t(0,0) - dot_product (p%x, t%t(1:,0))
    pt%x(1) = p%t * t%t(0,1) - dot_product (p%x, t%t(1:,1))
    pt%x(2) = p%t * t%t(0,2) - dot_product (p%x, t%t(1:,2))
    pt%x(3) = p%t * t%t(0,3) - dot_product (p%x, t%t(1:,3))
  end function momentum_tensor
end module omega_tensors
