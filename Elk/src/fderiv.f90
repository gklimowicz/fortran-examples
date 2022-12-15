
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: fderiv
! !INTERFACE:
subroutine fderiv(m,n,x,f,g)
! !INPUT/OUTPUT PARAMETERS:
!   m : order of derivative (in,integer)
!   n : number of points (in,integer)
!   x : abscissa array (in,real(n))
!   f : function array (in,real(n))
!   g : (anti-)derivative of f (out,real(n))
! !DESCRIPTION:
!   Given function $f$ defined on a set of points $x_i$ then if $m\ge 0$ this
!   routine computes the $m$th derivative of $f$ at each point. If $m=-1$ the
!   anti-derivative of $f$ given by
!   $$ g(x_i)=\int_{x_1}^{x_i} f(x)\,dx $$
!   is calculated. Both derivatives and integrals are computed by first fitting
!   the function to a clamped cubic spline.
!
! !REVISION HISTORY:
!   Created May 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m,n
real(8), intent(in) :: x(n),f(n)
real(8), intent(out) :: g(n)
! local variables
integer i
real(8) sm,dx
! automatic arrays
real(8) cf(3,n)
if (n.le.0) then
  write(*,*)
  write(*,'("Error(fderiv): invalid number of points : ",I8)') n
  write(*,*)
  stop
end if
! high accuracy integration/differentiation from spline interpolation
call spline(n,x,f,cf)
select case(m)
case(:-1)
  sm=0.d0
  g(1)=0.d0
  do i=1,n-1
    dx=x(i+1)-x(i)
    sm=sm+dx*(f(i) &
         +dx*(0.5d0*cf(1,i) &
         +dx*(0.3333333333333333333d0*cf(2,i) &
         +dx*0.25d0*cf(3,i))))
    g(i+1)=sm
  end do
case(1)
  g(:)=cf(1,:)
case(2)
  g(:)=2.d0*cf(2,:)
case(3)
  g(:)=6.d0*cf(3,:)
case(4:)
  g(:)=0.d0
end select
end subroutine
!EOC

