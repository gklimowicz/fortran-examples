
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine ztorfmt(nr,nri,zfmt,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: zfmt(*)
real(8), intent(out) :: rfmt(*)
! local variables
integer i
call ztorflmn(lmaxi,nri,lmmaxi,zfmt,rfmt)
i=lmmaxi*nri+1
call ztorflmn(lmaxo,nr-nri,lmmaxo,zfmt(i),rfmt(i))
return

contains

!BOP
! !ROUTINE: ztorflmn
! !INTERFACE:
pure subroutine ztorflmn(lmax,n,ld,zflm,rflm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   n    : number of functions to convert (in,integer)
!   ld   : leading dimension (in,integer)
!   zflm : coefficients of complex spherical harmonic expansion
!          (in,complex(ld,n))
!   rflm : coefficients of real spherical harmonic expansion (out,real(ld,n))
! !DESCRIPTION:
!   Converts a real function, $z_{lm}$, expanded in terms of complex spherical
!   harmonics into a real spherical harmonic expansion, $r_{lm}$:
!   $$ r_{lm}=\begin{cases}\frac{1}{\sqrt{2}}\Re(z_{lm}+(-1)^m z_{l-m}) & m>0 \\
!    \frac{1}{\sqrt{2}}\Im(-z_{lm}+(-1)^m z_{l-m}) & m<0 \\
!    \Re(z_{lm}) & m=0 \end{cases}\;. $$
!   See routine {\tt genrlm}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax,n,ld
complex(8), intent(in) :: zflm(ld,n)
real(8), intent(out) :: rflm(ld,n)
! local variables
integer l,m,lm1,lm2
! real constant 1/sqrt(2)
real(8), parameter :: c1=0.7071067811865475244d0
lm1=0
do l=0,lmax
  lm2=lm1+2*(l+1)
  do m=-l,-1
    lm1=lm1+1
    lm2=lm2-1
    if (mod(m,2).ne.0) then
      rflm(lm1,:)=-c1*(aimag(zflm(lm1,:))+aimag(zflm(lm2,:)))
    else
      rflm(lm1,:)=c1*(aimag(zflm(lm2,:))-aimag(zflm(lm1,:)))
    end if
  end do
  lm1=lm1+1
  lm2=lm2-1
  rflm(lm1,:)=dble(zflm(lm1,:))
  do m=1,l
    lm1=lm1+1
    lm2=lm2-1
    if (mod(m,2).ne.0) then
      rflm(lm1,:)=c1*(dble(zflm(lm1,:))-dble(zflm(lm2,:)))
    else
      rflm(lm1,:)=c1*(dble(zflm(lm1,:))+dble(zflm(lm2,:)))
    end if
  end do
end do
end subroutine
!EOC

end subroutine

