
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine rtozfmt(nr,nri,rfmt,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: rfmt(*)
complex(8), intent(out) :: zfmt(*)
! local variables
integer i
call rtozflmn(lmaxi,nri,lmmaxi,rfmt,zfmt)
i=lmmaxi*nri+1
call rtozflmn(lmaxo,nr-nri,lmmaxo,rfmt(i),zfmt(i))
return

contains

!BOP
! !ROUTINE: rtozflmn
! !INTERFACE:
pure subroutine rtozflmn(lmax,n,ld,rflm,zflm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   n    : number of functions to convert (in,integer)
!   ld   : leading dimension (in,integer)
!   rflm : coefficients of real spherical harmonic expansion (in,real(ld,n))
!   zflm : coefficients of complex spherical harmonic expansion
!          (out,complex(ld,n))
! !DESCRIPTION:
!   Converts a real function, $r_{lm}$, expanded in terms of real spherical
!   harmonics into a complex spherical harmonic expansion, $z_{lm}$:
!   $$ z_{lm}=\begin{cases} \frac{1}{\sqrt{2}}(r_{lm}+i(-1)^mr_{l-m}) & m>0 \\
!    \frac{1}{\sqrt{2}}((-1)^mr_{l-m}-ir_{lm}) & m<0 \\
!    r_{lm} & m=0 \end{cases}\;. $$
!   See routine {\tt genrlm}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax,n,ld
real(8), intent(in) :: rflm(ld,n)
complex(8), intent(out) :: zflm(ld,n)
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
      zflm(lm1,:)=c1*cmplx(-rflm(lm2,:),-rflm(lm1,:),8)
    else
      zflm(lm1,:)=c1*cmplx(rflm(lm2,:),-rflm(lm1,:),8)
    end if
  end do
  lm1=lm1+1
  lm2=lm2-1
  zflm(lm1,:)=rflm(lm1,:)
  do m=1,l
    lm1=lm1+1
    lm2=lm2-1
    if (mod(m,2).ne.0) then
      zflm(lm1,:)=c1*cmplx(rflm(lm1,:),-rflm(lm2,:),8)
    else
      zflm(lm1,:)=c1*cmplx(rflm(lm1,:),rflm(lm2,:),8)
    end if
  end do
end do
end subroutine
!EOC

end subroutine

