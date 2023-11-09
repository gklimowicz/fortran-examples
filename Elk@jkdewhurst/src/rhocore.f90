
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhocore
! !INTERFACE:
subroutine rhocore
! !USES:
use modmain
! !DESCRIPTION:
!   Adds the core density and magnetisation to the muffin-tin functions. Also
!   computes the amount of leakage of core charge from the muffin-tin spheres
!   into the interstitial.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Fixed core moment direction, October 2012 (M. Meinert)
!EOP
!BOC
implicit none
! local variables
integer ispn,idm,is,ias
integer nr,nri,iro,ir,i,i0,i1
real(8) v(ndmag),sm,t1
! external functions
real(8), external :: rfmtint
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
  sm=0.d0
! loop over spin channels
  do ispn=1,nspncr
! add the core density to the muffin-tin density
    i1=lmmaxi*(nri-1)+1
    rhomt(1:i1:lmmaxi,ias)=rhomt(1:i1:lmmaxi,ias)+rhocr(1:nri,ias,ispn)
    i0=i1+lmmaxi
    i1=lmmaxo*(nr-iro)+i0
    rhomt(i0:i1:lmmaxo,ias)=rhomt(i0:i1:lmmaxo,ias)+rhocr(iro:nr,ias,ispn)
! compute the core charge inside the muffin-tins
    t1=dot_product(wrmt(1:nr,is),rhocr(1:nr,ias,ispn))
    sm=sm+fourpi*y00*t1
  end do
! core leakage charge
  chgcrlk(ias)=chgcr(is)-sm
! add to the magnetisation in the case of a spin-polarised core
  if (spincore) then
! compute the moment in the muffin-tin
    do idm=1,ndmag
      v(idm)=rfmtint(nr,nri,wrmt(:,is),magmt(:,ias,idm))
    end do
! normalise
    if (ncmag) then
      t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
    else
      t1=abs(v(1))
    end if
    if (t1.gt.1.d-10) v(:)=v(:)/t1
! add the core magnetisation to the total
    i=1
    do ir=1,nri
      t1=rhocr(ir,ias,1)-rhocr(ir,ias,2)
      magmt(i,ias,:)=magmt(i,ias,:)+t1*v(:)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      t1=rhocr(ir,ias,1)-rhocr(ir,ias,2)
      magmt(i,ias,:)=magmt(i,ias,:)+t1*v(:)
      i=i+lmmaxo
    end do
  end if
end do
end subroutine
!EOC

