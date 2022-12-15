
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zpotclmt
! !INTERFACE:
pure subroutine zpotclmt(nr,nri,ld,rl,wpr,zrhomt,zvclmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr     : number of radial mesh points (in,integer)
!   nri    : number of points on inner part of muffin-tin (in,integer)
!   ld     : leading dimension (in,integer)
!   rl     : r^l on the radial mesh (in,real(ld,-lmaxo-1:lmaxo+2))
!   wpr    : weights for partial integration on radial mesh (in,real(4,nr))
!   zrhomt : muffin-tin charge density (in,complex(*))
!   zvclmt : muffin-tin Coulomb potential (out,complex(*))
! !DESCRIPTION:
!   Solves the Poisson equation for the charge density contained in an isolated
!   muffin-tin using the Green's function approach. In other words, the
!   spherical harmonic expansion of the Coulomb potential, $V_{lm}$, is obtained
!   from the density expansion, $\rho_{lm}$, by
!   $$ V_{lm}(r)=\frac{4\pi}{2l+1}\left(\frac{1}{r^{l+1}}\int_0^r\rho_{lm}(r')
!      {r'}^{l+2}dr'+r^l\int_r^R\frac{\rho_{lm}(r')}{{r'}^{l-1}}dr'\right) $$
!   where $R$ is the muffin-tin radius.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri,ld
real(8), intent(in) :: rl(ld,-lmaxo-1:lmaxo+2),wpr(4,nr)
complex(8), intent(in) :: zrhomt(*)
complex(8), intent(out) :: zvclmt(*)
! local variables
integer nro,iro,ir
integer l,l1,l2,l3
integer lm,npi,i0,i
real(8) t0
complex(8) z1
! automatic arrays
complex(8) f1(nr),f2(nr),f3(nr)
nro=nr-nri
iro=nri+1
npi=lmmaxi*nri
do l=0,lmaxi
  l1=l+2
  l2=-l+1
  l3=-l-1
  t0=fourpi/dble(2*l+1)
  do lm=l**2+1,(l+1)**2
    do ir=1,nri
      i=lm+lmmaxi*(ir-1)
      f1(ir)=rl(ir,l1)*zrhomt(i)
      f2(ir)=rl(ir,l2)*zrhomt(i)
    end do
    i0=lm+npi
    do ir=iro,nr
      i=i0+lmmaxo*(ir-iro)
      f1(ir)=rl(ir,l1)*zrhomt(i)
      f2(ir)=rl(ir,l2)*zrhomt(i)
    end do
    call splintwp(nr,wpr,f1,f3)
    call splintwp(nr,wpr,f2,f1)
    z1=f1(nr)
    do ir=1,nri
      i=lm+lmmaxi*(ir-1)
      zvclmt(i)=t0*(rl(ir,l3)*f3(ir)+rl(ir,l)*(z1-f1(ir)))
    end do
    do ir=iro,nr
      i=i0+lmmaxo*(ir-iro)
      zvclmt(i)=t0*(rl(ir,l3)*f3(ir)+rl(ir,l)*(z1-f1(ir)))
    end do
  end do
end do
do l=lmaxi+1,lmaxo
  l1=l+2
  l2=-l+1
  l3=-l-1
  t0=fourpi/dble(2*l+1)
  do lm=l**2+1,(l+1)**2
    i0=lm+npi
    do ir=iro,nr
      i=i0+lmmaxo*(ir-iro)
      f1(ir)=rl(ir,l1)*zrhomt(i)
      f2(ir)=rl(ir,l2)*zrhomt(i)
    end do
    call splintwp(nro,wpr(1,iro),f1(iro),f3(iro))
    call splintwp(nro,wpr(1,iro),f2(iro),f1(iro))
    z1=f1(nr)
    do ir=iro,nr
      i=i0+lmmaxo*(ir-iro)
      zvclmt(i)=t0*(rl(ir,l3)*f3(ir)+rl(ir,l)*(z1-f1(ir)))
    end do
  end do
end do
return

contains

pure subroutine splintwp(n,wp,f,g)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wp(*)
complex(8), intent(in) :: f(n)
complex(8), intent(out) :: g(n)
! local variables
integer i,j
complex(8) zsm
g(1)=0.d0
zsm=wp(5)*f(1)+wp(6)*f(2)+wp(7)*f(3)+wp(8)*f(4)
g(2)=zsm
do i=2,n-2
  j=i*4+1
  zsm=zsm+wp(j)*f(i-1)+wp(j+1)*f(i)+wp(j+2)*f(i+1)+wp(j+3)*f(i+2)
  g(i+1)=zsm
end do
j=(n-1)*4+1
g(n)=zsm+wp(j)*f(n-3)+wp(j+1)*f(n-2)+wp(j+2)*f(n-1)+wp(j+3)*f(n)
end subroutine

end subroutine
!EOC

