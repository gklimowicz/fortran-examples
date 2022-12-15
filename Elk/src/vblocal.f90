
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vblocal(vmt,vir,bmt)
use modmain
use modomp
implicit none
! arguments
real(8), intent(out) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(out) :: bmt(npcmtmax,natmtot,ndmag)
! local variables
integer idm,is,ias,nthd
integer nrc,nrci,npc
! automatic arrays
real(8) rfmt(npcmtmax)
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nrc,nrci) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
! convert muffin-tin Kohn-Sham potential to spherical coordinates
  call rfmtftoc(nrc,nrci,vsmt(:,ias),rfmt)
  call rbsht(nrc,nrci,rfmt,vmt(:,ias))
! multiply by radial integration weights
  call rfcmtwr(nrc,nrci,wrcmt(:,is),vmt(:,ias))
end do
!$OMP END PARALLEL DO
! multiply interstitial Kohn-Sham potential by characteristic function
vir(:)=vsir(:)*cfunir(:)
! repeat for the Kohn-Sham magnetic field
if (spinpol) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is,nrc,nrci,npc,idm) &
!$OMP NUM_THREADS(nthd)
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
    do idm=1,ndmag
      bmt(1:npc,ias,idm)=bsmt(1:npc,ias,idm)
      call rfcmtwr(nrc,nrci,wrcmt(:,is),bmt(:,ias,idm))
    end do
  end do
!$OMP END PARALLEL DO
end if
call freethd(nthd)
end subroutine

