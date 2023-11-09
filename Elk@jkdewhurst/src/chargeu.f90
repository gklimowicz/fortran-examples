
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine chargeu
use modmain
use modulr
use modtest
implicit none
! local variables
integer ifq,is,ias
integer nrc,nrci
real(8) t1
! automatic arrays
real(8) rfft(nqpt)
complex(8) zfft(nfqrz)
! external functions
real(8), external :: ddot
complex(8), external :: zfmtint
! calculate muffin-tin charges
chgmttot=0.d0
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  do ifq=1,nfqrz
    zfft(ifq)=zfmtint(nrc,nrci,wrcmt(:,is),rhoqmt(:,ias,ifq))
  end do
  chgmt(ias)=dble(zfft(1))
  chgmttot=chgmttot+chgmt(ias)
  call rzfftifc(3,ngridq,1,rfft,zfft)
  chgmtru(ias,:)=rfft(:)
end do
! calculate interstitial charge
t1=ddot(ngtc,rhoqir(:,1),2,cfrc,1)
chgir=t1*omega/dble(ngtc)
! total calculated charge
chgcalc=chgmttot+chgir
! write muffin-tin charges to file
call writetest(730,'ULR muffin-tin charges',nv=natmtot*nqpt,tol=5.d-5, &
 rva=chgmtru)
end subroutine

