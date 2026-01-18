
! Copyright (C) 2022 Leon Kerber, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tm3rtoz(l,k,p,r,ld,wkpr,zkpr)
use modmain
implicit none
! arguments
integer, intent(in) :: l,k,p,r,ld
real(8), intent(in) :: wkpr(-ld:ld)
complex(8), intent(out) :: zkpr(-ld:ld)
! local variables
integer g,t
real(8) a,b,t0,t1,t2
! external functions
real(8), external :: wigner3j,factn,factn2,factr
! old convention normalisation factors and phase factors
g=k+p+r
if (mod(g,2).eq.0) then
  t0=1.d0/wigner3j(k,p,r,0,0,0)
else
  t0=sqrt(factr(g+1,g-2*k)/(factn(g-2*p)*factn(g-2*r)))
  t0=t0*factn2(g-2*k)*factn2(g-2*p)*factn2(g-2*r)/factn2(g)
end if
t0=t0/sqrt(dble(2*r+1))
t0=t0*sqrt(factn(2*l-k)*factn(2*l+k+1))/factn(2*l)
t0=t0*sqrt(factn(2+p))
if (mod(k+p,2).ne.0) t0=-t0
! remove orthonormal convention normalisation factors
t0=t0/sqrt(dble(2*k+1))
t0=t0/sqrt(dble(2*p+1))
t0=t0/2.d0
do t=-r,r
  t1=t0*(wkpr(t)+wkpr(-t))
  t2=t0*(wkpr(t)-wkpr(-t))
  if (mod(t,2).eq.0) then
    a=t1
    b=t2
  else
    a=-t2
    b=-t1
  end if
  if ((k.eq.r).and.(p.eq.1)) then
    if (mod(k,2).eq.0) then
      b=-b
    else
      a=-a
    end if
  end if
  zkpr(t)=cmplx(a,b,8)
end do
end subroutine

