
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potuinit
use modmain
use modulr
use modrandom
use modomp
implicit none
! local variables
integer ifq,idm,is,ias
integer nrc,nrci,npc
real(8) cb,t1
! automatic arrays
real(8) rfmt1(npcmtmax),rfmt2(npcmtmax)
! set the Q=0 muffin-tin potential equal to that of the normal ground-state
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  call rfmtftoc(nrc,nrci,vsmt(:,ias),rfmt1)
  call rbsht(nrc,nrci,rfmt1,rfmt2)
  vsqmt(1:npc,ias,1)=rfmt2(1:npc)
end do
! zero the muffin-tin potential for non-zero Q
do ifq=2,nfqrz
  do ias=1,natmtot
    is=idxis(ias)
    vsqmt(1:npcmt(is),ias,ifq)=0.d0
  end do
end do
! repeat for the interstitial potential
vsqir(:,1)=vsir(:)
vsqir(:,2:nfqrz)=0.d0
if (.not.spinpol) return
! set the Q=0 muffin-tin magnetic field equal to that of the normal ground-state
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
    bsqmt(1:npc,ias,idm,1)=bsmt(1:npc,ias,idm)
  end do
end do
! zero the magnetic field for non-zero Q
do ifq=2,nfqrz
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      bsqmt(1:npcmt(is),ias,idm,ifq)=0.d0
    end do
  end do
end do
! repeat for the interstitial magnetic field
bsqir(:,:,1)=bsir(:,:)
bsqir(:,:,2:nfqrz)=0.d0
! coupling constant of the external field (g_e/4c)
cb=gfacte/(4.d0*solsc)
! initialise the external magnetic fields
t1=cb*rndbfcu
do ifq=1,nfqrz
  do idm=1,ndmag
    bfcq(idm,ifq)=t1*cmplx(randomu()-0.5d0,randomu()-0.5d0,8)
    do ias=1,natmtot
      bfcmtq(ias,idm,ifq)=t1*cmplx(randomu()-0.5d0,randomu()-0.5d0,8)
    end do
  end do
end do
bfcq(:,1)=dble(bfcq(:,1))
bfcmtq(:,:,1)=dble(bfcmtq(:,:,1))
end subroutine

