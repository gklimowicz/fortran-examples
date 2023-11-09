
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine plotulr(np,vpl,nf,zfmt,zfir,fp)
use modmain
use modulr
use modomp
implicit none
! arguments
integer, intent(in) :: np
real(8), intent(in) :: vpl(3,np)
integer, intent(in) :: nf
complex(8), intent(in) :: zfmt(npcmtmax,natmtot,nf,nfqrz)
complex(8), intent(in) :: zfir(ngtot,nf,nfqrz)
real(8), intent(out) :: fp(np,nf)
! local variables
integer iq,ifq0,ifq
integer jf,ip
real(8) sm,t1
complex(8) z1
! allocatable arrays
complex(8), allocatable :: fpq(:,:)
allocate(fpq(np,nfqrz))
! include or exclude the Q=0 component as required
if (tplotq0) then
  ifq0=1
else
  ifq0=2
end if
! loop over the number of functions
do jf=1,nf
! loop over real-complex FFT points
  do ifq=ifq0,nfqrz
! evaluate the complex function at all the plot points
    call zfplot(np,vpl,zfmt(:,:,jf,ifq),zfir(:,jf,ifq),fpq(:,ifq))
  end do
  do ip=1,np
    sm=0.d0
    do ifq=ifq0,nfqrz
      iq=iqrzf(ifq)
! multiply complex function by phase factor exp(iQ.r)
      t1=twopi*dot_product(vql(:,iq),vpl(:,ip))
      z1=cmplx(cos(t1),sin(t1),8)
      t1=dble(fpq(ip,ifq)*z1)
      if (ifq.gt.1) t1=t1*2.d0
      sm=sm+t1
    end do
    fp(ip,jf)=sm
  end do
end do
deallocate(fpq)
end subroutine

