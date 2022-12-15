
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeftm
use modmain
use moddftu
implicit none
! local variables
integer is,ia,ias,i
integer l,k,p,r,t
real(8) t0
! automatic arrays
real(8) wkpr(-lmmaxdm:lmmaxdm)
! open FTM.OUT
open(50,file='FTM.OUT',form='FORMATTED',action='WRITE')
do i=1,ntmfix
  is=itmfix(1,i)
  ia=itmfix(2,i)
  ias=idxas(ia,is)
  l=itmfix(3,i)
  k=itmfix(4,i)
  p=itmfix(5,i)
  r=itmfix(6,i)
  t=itmfix(7,i)
  write(50,*)
  write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
  write(50,'(" l = ",I2)') l
! scale factor for conventional normalisation
  t0=sqrt(dble((2*l+1)*nspinor))
  write(50,'(" k = ",I2,", p = ",I2,", r = ",I2,", t = ",I2)') k,p,r,t
! decompose density matrix in 3-index tensor moment components
  call dmtotm3(l,k,p,r,lmmaxdm,dmatmt(:,:,:,:,ias),wkpr)
  write(50,'(" tensor moment")')
  write(50,'("  current : ",G18.10)') t0*wkpr(t)
  write(50,'("  target  : ",G18.10)') wkprfix(i)
end do
close(50)
end subroutine

