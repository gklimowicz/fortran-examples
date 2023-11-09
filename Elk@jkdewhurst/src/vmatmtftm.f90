
! Copyright (C) 2014 L. Nordstrom, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vmatmtftm
use modmain
use moddftu
implicit none
! local variables
integer is,ia,ias,i
integer l,k,p,r,t
real(8) t0,t1
! automatic arrays
real(8) wkpr(-lmmaxdm:lmmaxdm)
complex(8) dm(lmmaxdm,2,lmmaxdm,2)
! loop over fixed tensor moment entries
do i=1,ntmfix
  is=itmfix(1,i)
  ia=itmfix(2,i)
  ias=idxas(ia,is)
  l=itmfix(3,i)
  k=itmfix(4,i)
  p=itmfix(5,i)
  r=itmfix(6,i)
  t=itmfix(7,i)
! decompose density matrix in 3-index tensor moment components
  call dmtotm3(l,k,p,r,lmmaxdm,dmatmt(:,:,:,:,ias),wkpr)
! scale factor for conventional normalisation
  t0=sqrt(dble((2*l+1)*2))
! take difference between current and target moment
  t1=wkpr(t)-wkprfix(i)/t0
  wkpr(:)=0.d0
  wkpr(t)=tauftm*t1
! compute new density matrix
  call tm3todm(l,k,p,r,lmmaxdm,wkpr,dm)
! add to global FTM potential matrix
  vmftm(:,:,:,:,ias)=vmftm(:,:,:,:,ias)+dm(:,:,:,:)
end do
! add to muffin-tin potential matrix (fix by Leon Kerber)
do i=1,ntmfix
  is=itmfix(1,i)
  ia=itmfix(2,i)
  ias=idxas(ia,is)
  vmatmt(:,:,:,:,ias)=vmatmt(:,:,:,:,ias)+vmftm(:,:,:,:,ias)
end do
end subroutine

