
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vblocalu
use modmain
use modulr
use modomp
implicit none
! local variables
integer ifq,idm,is,ias
integer nrc,nrci,npc,nthd
! automatic arrays
real(8) rfmt1(npcmtmax),rfmt2(npcmtmax)
! subtract the normal Kohn-Sham potential for Q=0
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  call rfmtftoc(nrc,nrci,vsmt(:,ias),rfmt1)
  call rbsht(nrc,nrci,rfmt1,rfmt2)
  vsqmt(1:npc,ias,1)=vsqmt(1:npc,ias,1)-rfmt2(1:npc)
end do
vsqir(1:ngtot,1)=vsqir(1:ngtot,1)-vsir(1:ngtot)
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ias,is,nrc,nrci) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ifq=1,nfqrz
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
! multiply the muffin-tin potential by the radial integration weights
    call zfcmtwr(nrc,nrci,wrcmt(:,is),vsqmt(:,ias,ifq))
  end do
end do
!$OMP END DO NOWAIT
!$OMP DO
do ifq=1,nfqrz
! multiply the interstitial potential by the characteristic function
  vsqir(1:ngtot,ifq)=vsqir(1:ngtot,ifq)*cfunir(1:ngtot)
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
if (.not.spinpol) return
! subtract the normal Kohn-Sham magnetic field for Q=0 in the muffin-tins
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
    bsqmt(1:npc,ias,idm,1)=bsqmt(1:npc,ias,idm,1)-bsmt(1:npc,ias,idm)
  end do
end do
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(idm,ias,is,nrc,nrci) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ifq=1,nfqrz
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      nrc=nrcmt(is)
      nrci=nrcmti(is)
! multiply the muffin-tin field by the radial integration weights
      call zfcmtwr(nrc,nrci,wrcmt(:,is),bsqmt(:,ias,idm,ifq))
    end do
  end do
end do
!$OMP END DO NOWAIT
!$OMP DO
do ifq=1,nfqrz
  do idm=1,ndmag
! multiply interstitial field by the characteristic function
    bsqir(1:ngtot,idm,ifq)=bsqir(1:ngtot,idm,ifq)*cfunir(1:ngtot)
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! subtract the normal Kohn-Sham magnetic field for Q=0 in the interstitial
! (this is already multiplied by the characteristic function)
do idm=1,ndmag
  bsqir(1:ngtot,idm,1)=bsqir(1:ngtot,idm,1)-bsir(1:ngtot,idm)
end do
end subroutine

