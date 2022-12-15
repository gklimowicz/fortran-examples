
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine momentu
use modmain
use modulr
use modtest
implicit none
! local variables
integer ifq,idm,is,ias
integer nrc,nrci
real(8) t1
! automatic arrays
real(8) rfft(nqpt)
complex(8) zfft(nfqrz)
! external functions
real(8), external :: ddot
complex(8), external :: zfmtint
if (.not.spinpol) return
! calculate muffin-tin moments
mommttot(:)=0.d0
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    do ifq=1,nfqrz
      zfft(ifq)=zfmtint(nrc,nrci,wrcmt(:,is),magqmt(:,ias,idm,ifq))
    end do
    mommt(idm,ias)=dble(zfft(1))
    mommttot(idm)=mommttot(idm)+mommt(idm,ias)
    call rzfftifc(3,ngridq,1,rfft,zfft)
    mommtru(idm,ias,:)=rfft(:)
  end do
end do
! calculate interstitial moment
do idm=1,ndmag
  t1=ddot(ngtc,magqir(:,idm,1),2,cfrc,1)
  momir(idm)=t1*omega/dble(ngtc)
end do
momtot(:)=mommttot(:)+momir(:)
! total moment magnitude
if (ncmag) then
  momtotm=sqrt(momtot(1)**2+momtot(2)**2+momtot(3)**2)
else
  momtotm=abs(momtot(1))
end if
! write the muffin-tin moments to test file
call writetest(770,'ULR muffin-tin moments',nv=ndmag*natmtot*nqpt,tol=1.d-2, &
 rva=mommtru)
end subroutine

