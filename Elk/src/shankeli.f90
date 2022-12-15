
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: shankeli
! !INTERFACE:
subroutine shankeli(lmax,x,hl)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum order of Hankel function (in,integer)
!   x    : real argument (in,real)
!   hl   : array of returned values (out,real(0:lmax))
! !DESCRIPTION:
!   Computes the spherical Hankel function of the first kind with imaginary
!   argument, $\tilde{h}_l(x)=i^lh_l(ix)$, for real $x$ and
!   $l=0\ldots l_{\rm max}$. The recurrence relation
!   $$ \tilde{h}_{l+1}(x)=\frac{2l+1}{x}\tilde{h}_l(x)+\tilde{h}_{l-1}(x) $$
!   is used upwards. The starting values there are
!   $\tilde{h}_0(x)=-e^{-x}/x$ and $\tilde{h}_1(x)=\tilde{h}_0(x)(1+1/x)$.
!   For $x\ll 1$ we use the asymptotic form
!   $$ \tilde{h}_l(x)\approx\frac{-(2l-1)!!}{(-x)^{l+1}}. $$
!
! !REVISION HISTORY:
!   Created April 2008 from sbessel routine (Lars Nordstrom)
!   Changed name, September 2021 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
real(8), intent(in) :: x
real(8), intent(out) :: hl(0:lmax)
! local variables
integer l
real(8) xi,h0,h1,t1
if ((lmax.lt.0).or.(lmax.gt.50)) then
  write(*,*)
  write(*,'("Error(shankeli): lmax out of range : ",I8)') lmax
  write(*,*)
  stop
end if
if ((x.le.0.d0).or.(x.gt.1.d8)) then
  write(*,*)
  write(*,'("Error(shankeli): x out of range : ",G18.10)') x
  write(*,*)
  stop
end if
xi=1.d0/x
hl(0)=-xi*exp(-x)
if (lmax.eq.0) return
! treat x << 1
if (x.lt.1.d-8) then
  t1=-xi
  do l=1,lmax
    t1=t1*xi*dble(2*l-1)
    hl(l)=t1
  end do
  return
end if
! recurse up
hl(1)=hl(0)*(1.d0+xi)
if (lmax.eq.1) return
h0=hl(0)
h1=hl(1)
do l=2,lmax
  t1=(2*l-1)*h1*xi+h0
  h0=h1
  h1=t1
  hl(l)=h1
end do
end subroutine
!EOC

