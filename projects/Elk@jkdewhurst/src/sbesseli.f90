
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: sbesseli
! !INTERFACE:
subroutine sbesseli(lmax,x,jl)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum order of Bessel function (in,integer)
!   x    : real argument (in,real)
!   jl   : array of returned values (out,real(0:lmax))
! !DESCRIPTION:
!   Computes spherical Bessel functions with imaginary argument,
!   $\tilde{j}_l(x)\equiv i^lj_l(ix)$, for real $x$ and
!   $l=0\ldots l_{\rm max}$. The recurrence relation
!   $$ \tilde{j}_{l+1}(x)=\frac{2l+1}{x}\tilde{j}_l(x)+\tilde{j}_{l-1}(x) $$
!   is used either downwards for $x<2\,l_{\rm max}$ or upwards for
!   $x\ge 2\,l_{\rm max}$. The starting values are $\tilde{j}_0(x)=\sinh(x)/x$
!   and $\tilde{j}_1(x)=(\tilde{j}_0(x)-\cosh(x))/x$. The asymptotic form
!   $$ \tilde{j}_l(x)\approx\frac{(-x)^l}{(2l+1)!!} $$
!   is used for $x\ll 1$.
!
! !REVISION HISTORY:
!   Created April 2008 from sbessel routine (Lars Nordstrom)
!   Fixed accuracy issue and changed name, September 2021 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
real(8), intent(in) :: x
real(8), intent(out) :: jl(0:lmax)
! local variables
integer l,lst
real(8), parameter :: rsc=1.d150,rsci=1.d0/rsc
real(8) xi,j0,j1,t1
if ((lmax.lt.0).or.(lmax.gt.20)) then
  write(*,*)
  write(*,'("Error(sbesseli): lmax out of range : ",I8)') lmax
  write(*,*)
  stop
end if
if ((x.lt.0.d0).or.(x.gt.1.d8)) then
  write(*,*)
  write(*,'("Error(sbesseli): x out of range : ",G18.10)') x
  write(*,*)
  stop
end if
! treat x << 1
if (x.lt.1.d-8) then
  jl(0)=1.d0
  t1=1.d0
  do l=1,lmax
    t1=-t1*x/dble(2*l+1)
    jl(l)=t1
  end do
  return
end if
if (lmax.eq.0) then
  jl(0)=sinh(x)/x
  return
end if
xi=1.d0/x
if (x.lt.2*lmax) then
! for x < 2*lmax recurse down
  j1=1.d0
  j0=0.d0
! starting value for l above lmax
  lst=lmax+lmax/2+12
  do l=lst,lmax+1,-1
    t1=j0-(2*l+1)*j1*xi
    j0=j1
    j1=t1
! check for overflow
    if (abs(j1).gt.rsc) then
! rescale
      t1=t1*rsci
      j0=j0*rsci
      j1=j1*rsci
    end if
  end do
  do l=lmax,0,-1
    t1=j0-(2*l+1)*j1*xi
    j0=j1
    j1=t1
! check for overflow
    if (abs(j1).gt.rsc) then
! rescale
      t1=t1*rsci
      j0=j0*rsci
      j1=j1*rsci
      jl(l+1:lmax)=jl(l+1:lmax)*rsci
    end if
    jl(l)=j0
  end do
! rescaling constant
  t1=sinh(x)/(x*j0)
  jl(:)=t1*jl(:)
else
! for x >= 2*lmax recurse up
  jl(0)=sinh(x)*xi
  jl(1)=(jl(0)-cosh(x))*xi
  if (lmax.eq.1) return
  j0=jl(0)
  j1=jl(1)
  do l=2,lmax
    t1=(2*l-1)*j1*xi+j0
    j0=j1
    j1=t1
    jl(l)=j1
  end do
end if
end subroutine
!EOC

