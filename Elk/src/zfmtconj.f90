
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine zfmtconj(nr,nri,np,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri,np
complex(8), intent(inout) :: zfmt(np)
! local variables
integer i
! automatic arrays
complex(8) zfmt1(np)
zfmt1(:)=zfmt(:)
call zflmnconj(lmaxi,nri,lmmaxi,zfmt1,zfmt)
i=lmmaxi*nri+1
call zflmnconj(lmaxo,nr-nri,lmmaxo,zfmt1(i),zfmt(i))
return

contains

!BOP
! !ROUTINE: zflmnconj
! !INTERFACE:
pure subroutine zflmnconj(lmax,n,ld,zflm1,zflm2)
! !INPUT/OUTPUT PARAMETERS:
!   lmax  : maximum angular momentum (in,integer)
!   n     : number of functions to conjugate (in,integer)
!   ld    : leading dimension (in,integer)
!   zflm1 : coefficients of input complex spherical harmonic expansion
!           (in,complex((lmax+1)**2)))
!   zflm2 : coefficients of output complex spherical harmonic expansion
!           (out,complex((lmax+1)**2)))
! !DESCRIPTION:
!   Returns the complex conjugate of a function expanded in spherical harmonics.
!   In other words, given the input function coefficients $z_{lm}$, the routine
!   returns  $z'_{lm}=(-1)^m z^*_{l-m}$ so that
!   $$ \sum_{lm}z'_{lm}Y_{lm}(\theta,\phi)=\left(\sum_{lm}z_{lm}Y_{lm}
!    (\theta,\phi)\right)^* $$
!   for all $(\theta,\phi)$.
!
! !REVISION HISTORY:
!   Created April 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax,n,ld
complex(8), intent(in) :: zflm1(ld,n)
complex(8), intent(out) :: zflm2(ld,n)
! local variables
integer l,m,lm1,lm2
do l=0,lmax
  lm1=l**2
  lm2=(l+1)**2+1
  do m=-l,-1
    lm1=lm1+1
    lm2=lm2-1
    if (mod(m,2).eq.0) then
      zflm2(lm1,:)=conjg(zflm1(lm2,:))
      zflm2(lm2,:)=conjg(zflm1(lm1,:))
    else
      zflm2(lm1,:)=-conjg(zflm1(lm2,:))
      zflm2(lm2,:)=-conjg(zflm1(lm1,:))
    end if
  end do
! m=0 case
  lm1=lm1+1
  zflm2(lm1,:)=conjg(zflm1(lm1,:))
end do
end subroutine
!EOC

end subroutine

