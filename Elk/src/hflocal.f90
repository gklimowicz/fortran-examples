
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hflocal(vmt,vir,bmt,bir)
use modmain
implicit none
! arguments
real(8), intent(out) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(out) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
! local variables
integer idm,is,ias
integer np,nrc,nrci
! automatic arrays
real(8) rfmt1(npmtmax),rfmt2(npcmtmax)
! compute the Coulomb potential
call potcoul
! convert to spherical coordinates and store in output arrays
if (hybrid) then
! hybrid functional case
  call potxc(.true.,xctype,rhomt,rhoir,magmt,magir,taumt,tauir,exmt,exir,ecmt, &
   ecir,vxcmt,vxcir,bxcmt,bxcir,wxcmt,wxcir)
  do ias=1,natmtot
    is=idxis(ias)
    np=npmt(is)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    rfmt1(1:np)=vclmt(1:np,ias)+vxcmt(1:np,ias)
    call rfmtftoc(nrc,nrci,rfmt1,rfmt2)
    call rbsht(nrc,nrci,rfmt2,vmt(:,ias))
    call rfcmtwr(nrc,nrci,wrcmt(:,is),vmt(:,ias))
  end do
  vir(:)=(vclir(:)+vxcir(:))*cfunir(:)
  if (spinpol) then
    do idm=1,ndmag
      do ias=1,natmtot
        is=idxis(ias)
        nrc=nrcmt(is)
        nrci=nrcmti(is)
        call rfmtftoc(nrc,nrci,bxcmt(:,ias,idm),rfmt1)
        call rbsht(nrc,nrci,rfmt1,bmt(:,ias,idm))
        call rfcmtwr(nrc,nrci,wrcmt(:,is),bmt(:,ias,idm))
      end do
      bir(:,idm)=bxcir(:,idm)*cfunir(:)
    end do
  end if
else
! normal Hartree-Fock case
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    call rfmtftoc(nrc,nrci,vclmt(:,ias),rfmt1)
    call rbsht(nrc,nrci,rfmt1,vmt(:,ias))
    call rfcmtwr(nrc,nrci,wrcmt(:,is),vmt(:,ias))
  end do
  vir(:)=vclir(:)*cfunir(:)
end if
end subroutine

