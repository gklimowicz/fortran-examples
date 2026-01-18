
! Copyright (C) 2022 Leon Kerber, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetm3td
use modmain
use modtddft
use moddftu
implicit none
! local variables
integer is,ia,ias,idu
integer l,p,k,r,t,i
real(8) t0
character(256) fname
! automatic arrays
real(8) wkpr(-lmmaxdm:lmmaxdm)
! loop over DFT+U entries
i=0
do idu=1,ndftu
  is=isldu(1,idu)
  l=isldu(2,idu)
! scale factor for conventional normalisation
  t0=sqrt(dble((2*l+1)*2))
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do k=0,2*l
      do p=0,1
        do r=abs(k-p),k+p
          i=i+1
! decompose density matrix in 3-index tensor moment components
          call dmtotm3(l,k,p,r,lmmaxdm,dmatmt(:,:,:,:,ias),wkpr)
! construct the filename
          write(fname,'("TMTD_S",I2.2,"_A",I3.3,"_L",I1,"_K",I1,"_P",I1,&
           &"_R",I1,".OUT")') is,ia,l,k,p,r
          if (itimes.le.1) then
            open(50,file=trim(fname),form='FORMATTED')
          else
            open(50,file=trim(fname),form='FORMATTED',position='APPEND')
          end if
          write(50,'(20G18.10)') times(itimes),(t0*(wkpr(t)-wkpr0(t,i)),t=-r,r)
          close(50)
        end do
      end do
    end do
  end do
end do
end subroutine

