
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: genylmv
! !INTERFACE:
pure subroutine genylmv(lmax,v,ylm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   v    : input vector (in,real(3))
!   ylm  : array of spherical harmonics (out,complex((lmax+1)**2))
! !DESCRIPTION:
!   Generates a sequence of spherical harmonics, including the Condon-Shortley
!   phase, evaluated at angles $(\theta,\phi)$ for $0<l<l_{\rm max}$. The values
!   are returned in a packed array {\tt ylm} indexed with $j=l(l+1)+m+1$. This
!   routine is numerically stable and accurate to near machine precision for
!   $l\le 50$.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!   Improved stability, December 2005 (JKD)
!   Changed algorithm, June 2019 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
real(8), intent(in) :: v(3)
complex(8), intent(out) :: ylm(*)
! local variables
integer l,m,lm1,lm2,lm3,lm4
real(8), parameter :: eps=1.d-14
real(8) r,st,ct,sp,cp
real(8) t1,t2,t3,t4
complex(8) z1
ylm(1)=0.28209479177387814347d0
if (lmax.eq.0) return
r=sqrt(v(1)**2+v(2)**2+v(3)**2)
if (r.gt.eps) then
  t1=v(3)/r
  if (t1.ge.1.d0) then
    st=0.d0
    ct=1.d0
  else if (t1.le.-1.d0) then
    st=0.d0
    ct=-1.d0
  else
    st=sqrt(1.d0-t1**2)
    ct=t1
  end if
  if ((abs(v(1)).gt.eps).or.(abs(v(2)).gt.eps)) then
    t1=1.d0/sqrt(v(1)**2+v(2)**2)
    sp=t1*v(2)
    cp=t1*v(1)
  else
    sp=0.d0
    cp=1.d0
  end if
else
  st=0.d0
  ct=1.d0
  sp=0.d0
  cp=1.d0
end if
z1=cmplx(cp,sp,8)
ylm(3)=0.48860251190291992159d0*ct
ylm(4)=-0.34549414947133547927d0*st*z1
ylm(2)=-conjg(ylm(4))
do l=2,lmax
  lm1=(l+1)**2
  lm2=l**2
  lm3=(l-1)**2+1
  lm4=lm2+1
  ylm(lm1)=-st*sqrt(dble(2*l+1)/dble(2*l))*z1*ylm(lm2)
  if (mod(l,2).eq.0) then
    ylm(lm4)=conjg(ylm(lm1))
  else
    ylm(lm4)=-conjg(ylm(lm1))
  end if
  lm1=lm1-1
  ylm(lm1)=ct*sqrt(dble(2*l+1))*ylm(lm2)
  lm4=lm4+1
  if (mod(l-1,2).eq.0) then
    ylm(lm4)=conjg(ylm(lm1))
  else
    ylm(lm4)=-conjg(ylm(lm1))
  end if
  t1=ct*sqrt(dble((2*l-1)*(2*l+1)))
  t2=sqrt(dble((2*l+1))/dble(2*l-3))
  do m=l-2,1,-1
    lm1=lm1-1; lm2=lm2-1; lm3=lm3-1; lm4=lm4+1
    t3=1.d0/sqrt(dble((l-m)*(l+m)))
    t4=t2*sqrt(dble((l-m-1)*(l+m-1)))
    ylm(lm1)=t3*(t1*ylm(lm2)-t4*ylm(lm3))
    if (mod(m,2).eq.0) then
      ylm(lm4)=conjg(ylm(lm1))
    else
      ylm(lm4)=-conjg(ylm(lm1))
    end if
  end do
  lm1=lm1-1; lm2=lm2-1; lm3=lm3-1
  t3=1.d0/dble(l)
  t4=t2*dble(l-1)
  ylm(lm1)=t3*(t1*ylm(lm2)-t4*ylm(lm3))
end do
end subroutine
!EOC

