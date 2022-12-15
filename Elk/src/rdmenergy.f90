
! Copyright (C) 2005-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmenergy
! !INTERFACE:
subroutine rdmenergy
! !USES:
use modmain
use modrdm
use modtest
! !DESCRIPTION:
!   Calculates RDMFT total energy (free energy for finite temperatures).
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!   Updated for free energy 2009 (Baldsiefen)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,is,ias
integer nr,nri,ir,i
real(8) wo
complex(8) z1
! automatic arrays
real(8) rfmt(npmtmax)
! allocatable arrays
complex(8), allocatable :: evecsv(:,:)
! external functions
real(8), external :: rfmtinp
complex(8), external :: zdotc
! Coulomb energy from core states
engyvcl=0.d0
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  rfmt(1:npmt(is))=0.d0
  i=1
  if (spincore) then
! spin-polarised core
    do ir=1,nri
      rfmt(i)=rhocr(ir,ias,1)+rhocr(ir,ias,2)
      i=i+lmmaxi
    end do
    do ir=nri+1,nr
      rfmt(i)=rhocr(ir,ias,1)+rhocr(ir,ias,2)
      i=i+lmmaxo
    end do
  else
! spin-unpolarised core
    do ir=1,nri
      rfmt(i)=rhocr(ir,ias,1)
      i=i+lmmaxi
    end do
    do ir=nri+1,nr
      rfmt(i)=rhocr(ir,ias,1)
      i=i+lmmaxo
    end do
  end if
  engyvcl=engyvcl+rfmtinp(nr,nri,wrmt(:,is),rfmt,vclmt(:,ias))
end do
engykn=engykncr
allocate(evecsv(nstsv,nstsv))
do ik=1,nkpt
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
  do ist=1,nstsv
    wo=wkpt(ik)*occsv(ist,ik)
! Coulomb energy from valence states
    engyvcl=engyvcl+wo*dble(vclmat(ist,ist,ik))
! kinetic energy from valence states
    z1=zdotc(nstsv,evecsv(:,ist),1,dkdc(:,ist,ik),1)
    engykn=engykn+wo*dble(z1)
  end do
end do
deallocate(evecsv)
! Madelung term
engymad=0.d0
do ias=1,natmtot
  is=idxis(ias)
  engymad=engymad+0.5d0*spzn(is)*(vclmt(1,ias)-vcln(1,is))*y00
end do
! exchange-correlation energy
call rdmengyxc
! total energy
engytot=0.5d0*engyvcl+engymad+engykn+engyx
if (rdmtemp.gt.0.d0) then
  call rdmentropy
  engytot=engytot-rdmtemp*rdmentrpy
end if
! write total energy to test file
call writetest(300,'RDMFT total energy',tol=1.d-6,rv=engytot)
end subroutine
!EOC

