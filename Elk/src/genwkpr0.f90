
! Copyright (C) 2022 Leon Kerber, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genwkpr0
use modmain
use modtddft
use moddftu
implicit none
integer is,ia,ias,idu
integer n,l,p,k,r,i
! determine the total number of tensor moments
n=0
do idu=1,ndftu
  is=isldu(1,idu)
  l=isldu(2,idu)
  do ia=1,natoms(is)
    do k=0,2*l
      do p=0,1
        do r=abs(k-p),k+p
          n=n+1
        end do
      end do
    end do
  end do
end do
! allocate the t=0 tensor moment global array
if (allocated(wkpr0)) deallocate(wkpr0)
allocate(wkpr0(-lmmaxdm:lmmaxdm,n))
i=0
do idu=1,ndftu
  is=isldu(1,idu)
  l=isldu(2,idu)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do k=0,2*l
      do p=0,1
        do r=abs(k-p),k+p
          i=i+1
! decompose density matrix in 3-index tensor moment components
          call dmtotm3(l,k,p,r,lmmaxdm,dmatmt(:,:,:,:,ias),wkpr0(:,i))
        end do
      end do
    end do
  end do
end do
end subroutine

