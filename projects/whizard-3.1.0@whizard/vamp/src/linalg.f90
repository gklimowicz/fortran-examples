! linalg.f90 --
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
module linalg
  use kinds
  use utils
  implicit none
  private
  public :: lu_decompose
  public :: determinant
  public :: diagonalize_real_symmetric
  private :: jacobi_rotation
  public :: unit, diag
contains
  pure subroutine lu_decompose (a, pivots, eps, l, u)
    real(kind=default), dimension(:,:), intent(inout) :: a
    integer, dimension(:), intent(out), optional :: pivots
    real(kind=default), intent(out), optional :: eps
    real(kind=default), dimension(:,:), intent(out), optional :: l, u
    real(kind=default), dimension(size(a,dim=1)) :: vv
    integer, dimension(size(a,dim=1)) :: p
    integer :: j, pivot
    if (present (eps)) then
       eps = 1.0
    end if
    vv = maxval (abs (a), dim=2)
    if (any (vv == 0.0)) then
       a = 0.0
       if (present (pivots)) then
          pivots = 0
       end if
       if (present (eps)) then
          eps = 0
       end if
       return
    end if
    vv = 1.0 / vv
    do j = 1, size (a, dim=1)
       pivot = j - 1 + sum (maxloc (vv(j:) * abs (a(j:,j))))
       if (j /= pivot) then
          call swap (a(pivot,:), a(j,:))
          if (present (eps)) then
             eps = - eps
          end if
          vv(pivot) = vv(j)
       end if
       p(j) = pivot
       if (a(j,j) == 0.0) then
          a(j,j) = tiny (a(j,j))
       end if
       a(j+1:,j) = a(j+1:,j) / a(j,j)
       a(j+1:,j+1:) &
            = a(j+1:,j+1:) - outer_product (a(j+1:,j), a(j,j+1:))
    end do
    if (present (pivots)) then
       pivots = p
    end if
    if (present (l)) then
       do j = 1, size (a, dim=1)
          l(1:j-1,j) = 0.0
          l(j,j) = 1.0
          l(j+1:,j) = a(j+1:,j)
       end do
       do j = size (a, dim=1), 1, -1
          call swap (l(j,:), l(p(j),:))
       end do
    end if
    if (present (u)) then
       do j = 1, size (a, dim=1)
          u(1:j,j) = a(1:j,j)
          u(j+1:,j) = 0.0
       end do
    end if
  end subroutine lu_decompose
  pure subroutine determinant (a, det)
    real(kind=default), dimension(:,:), intent(in) :: a
    real(kind=default), intent(out) :: det
    real(kind=default), dimension(size(a,dim=1),size(a,dim=2)) :: lu
    integer :: i
    lu = a
    call lu_decompose (lu, eps = det)
    do i = 1, size (a, dim = 1)
       det = det * lu(i,i)
    end do
  end subroutine determinant
  pure subroutine diagonalize_real_symmetric (a, eval, evec, num_rot)
    real(kind=default), dimension(:,:), intent(in) :: a
    real(kind=default), dimension(:), intent(out) :: eval
    real(kind=default), dimension(:,:), intent(out) :: evec
    integer, intent(out), optional :: num_rot
    real(kind=default), dimension(size(a,dim=1),size(a,dim=2)) :: aa
    real(kind=default) :: off_diagonal_norm, threshold, &
         c, g, h, s, t, tau, cot_2phi
    logical, dimension(size(eval),size(eval)) :: upper_triangle
    integer, dimension(size(eval)) :: one_to_ndim
    integer :: p, q, ndim, j, sweep
    integer, parameter :: MAX_SWEEPS = 50
    ndim = size (eval)
    one_to_ndim = (/ (j, j=1,ndim) /)
    upper_triangle = &
         spread (one_to_ndim, dim=1, ncopies=ndim) &
           > spread (one_to_ndim, dim=2, ncopies=ndim)
    aa = a
    call unit (evec)
    if (present (num_rot)) then
       num_rot = 0
    end if
    sweeps: do sweep = 1, MAX_SWEEPS
       off_diagonal_norm = sum (abs (aa), mask=upper_triangle)
       if (off_diagonal_norm == 0.0) then
          eval = diag (aa)
          return
       end if
       if (sweep < 4) then
          threshold = 0.2 * off_diagonal_norm / ndim**2
       else
          threshold = 0.0
       end if
       do p = 1, ndim - 1
          do q = p + 1, ndim
             g = 100 * abs (aa (p,q))
             if ((sweep > 4) &
                  .and. (g <= min (spacing (aa(p,p)), spacing (aa(q,q))))) then
                aa(p,q) = 0.0
             else if (abs (aa(p,q)) > threshold) then
                h = aa(q,q) - aa(p,p)
                if (g <= spacing (h)) then
                   t = aa(p,q) / h
                else
                   cot_2phi = 0.5 * h / aa(p,q)
                   t = sign (1.0_default, cot_2phi) &
                         / (abs (cot_2phi) + sqrt (1.0 + cot_2phi**2))
                end if
                c = 1.0 / sqrt (1.0 + t**2)
                s = t * c
                tau = s / (1.0 + c)
                aa(p,p) = aa(p,p) - t * aa(p,q)
                aa(q,q) = aa(q,q) + t * aa(p,q)
                aa(p,q) = 0.0
                call jacobi_rotation (s, tau, aa(1:p-1,p), aa(1:p-1,q))
                call jacobi_rotation (s, tau, aa(p,p+1:q-1), aa(p+1:q-1,q))
                call jacobi_rotation (s, tau, aa(p,q+1:ndim), aa(q,q+1:ndim))
                call jacobi_rotation (s, tau, evec(:,p), evec(:,q))
                if (present (num_rot)) then
                   num_rot = num_rot + 1
                end if
             end if
          end do
       end do
    end do sweeps
    if (present (num_rot)) then
       num_rot = -1
    end if
  !!! print *, "linalg::diagonalize_real_symmetric: exceeded sweep count"
  end subroutine diagonalize_real_symmetric
  pure subroutine jacobi_rotation (s, tau, vp, vq)
    real(kind=default), intent(in) :: s, tau
    real(kind=default), dimension(:), intent(inout) :: vp, vq
    real(kind=default), dimension(size(vp)) :: vp_tmp
    vp_tmp = vp
    vp = vp - s * (vq     + tau * vp)
    vq = vq + s * (vp_tmp - tau * vq)
  end subroutine jacobi_rotation
  pure subroutine unit (u)
    real(kind=default), dimension(:,:), intent(out) :: u
    integer :: i
    u = 0.0
    do i = 1, min (size (u, dim = 1), size (u, dim = 2))
       u(i,i) = 1.0
    end do
  end subroutine unit
  pure function diag (a) result (d)
    real(kind=default), dimension(:,:), intent(in) :: a
    real(kind=default), dimension(min(size(a,dim=1),size(a,dim=2))) :: d
    integer :: i
    do i = 1, min (size (a, dim = 1), size (a, dim = 2))
       d(i) = a(i,i)
    end do
  end function diag
end module linalg
