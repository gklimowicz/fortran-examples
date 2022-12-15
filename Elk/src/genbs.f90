
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine genbs
use modmain
use modomp
implicit none
! local variables
integer idm,is,ia,ias
integer nrc,nrci,npc,nthd
real(8) cb,t1
! automatic arrays
real(8) rfmt(npcmtmax)
if (.not.spinpol) return
! coupling constant of the external field (g_e/4c)
cb=gfacte/(4.d0*solsc)
!------------------------------------!
!     muffin-tin Kohn-Sham field     !
!------------------------------------!
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,ia,nrc) &
!$OMP PRIVATE(nrci,npc,idm,t1) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
! exchange-correlation magnetic field in spherical coordinates
  do idm=1,ndmag
    call rfmtftoc(nrc,nrci,bxcmt(:,ias,idm),rfmt)
    call rbsht(nrc,nrci,rfmt,bsmt(:,ias,idm))
  end do
! add the external magnetic field
  t1=cb*(bfcmt(3,ia,is)+bfieldc(3))
  bsmt(1:npc,ias,ndmag)=bsmt(1:npc,ias,ndmag)+t1
  if (ncmag) then
    do idm=1,2
      t1=cb*(bfcmt(idm,ia,is)+bfieldc(idm))
      bsmt(1:npc,ias,idm)=bsmt(1:npc,ias,idm)+t1
    end do
  end if
end do
!$OMP END PARALLEL DO
call freethd(nthd)
!-----------------------------------------------!
!     interstitial Kohn-Sham magnetic field     !
!-----------------------------------------------!
call holdthd(ndmag,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(t1) &
!$OMP NUM_THREADS(nthd)
do idm=1,ndmag
  if (ncmag) then
    t1=cb*bfieldc(idm)
  else
    t1=cb*bfieldc(3)
  end if
! multiply by characteristic function
  bsir(1:ngtot,idm)=(bxcir(1:ngtot,idm)+t1)*cfunir(1:ngtot)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! add the magnetic dipole field if required
if (tbdip) call bdipole
end subroutine

