
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rhomagsh
! !INTERFACE:
subroutine rhomagsh
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Converts the muffin-tin density and magnetisation from spherical coordinates
!   to a spherical harmonic expansion. See {\tt rhomagk}.
!
! !REVISION HISTORY:
!   Created January 2009 (JKD)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ias
integer npc,nthd
! automatic arrays
real(8) rfmt(npcmtmax)
call holdthd(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,npc,idm) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  npc=npcmt(is)
! convert the density to spherical harmonics
  rfmt(1:npc)=rhomt(1:npc,ias)
  call rfsht(nrcmt(is),nrcmti(is),rfmt,rhomt(:,ias))
end do
!$OMP END DO NOWAIT
do idm=1,ndmag
!$OMP DO
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
! convert the magnetisation to spherical harmonics
    rfmt(1:npc)=magmt(1:npc,ias,idm)
    call rfsht(nrcmt(is),nrcmti(is),rfmt,magmt(:,ias,idm))
  end do
!$OMP END DO
end do
!$OMP END PARALLEL
call freethd(nthd)
end subroutine
!EOC

