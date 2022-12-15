
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writestate
! !INTERFACE:
subroutine writestate
! !USES:
use modmain
use moddftu
! !DESCRIPTION:
!   Writes the charge density, potentials and other relevant variables to the
!   file {\tt STATE.OUT}. Note to developers: changes to the way the variables
!   are written should be mirrored in {\tt readstate}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ias
! allocatable arrays
real(8), allocatable :: rfmt(:,:,:),rvfmt(:,:,:,:),rvfcmt(:,:,:,:)
open(100,file='STATE'//trim(filext),form='UNFORMATTED',action='WRITE')
write(100) version
write(100) spinpol
write(100) nspecies
write(100) lmmaxo
write(100) nrmtmax
write(100) nrcmtmax
do is=1,nspecies
  write(100) natoms(is)
  write(100) nrmt(is)
  write(100) rsp(1:nrmt(is),is)
  write(100) nrcmt(is)
  write(100) rcmt(1:nrcmt(is),is)
end do
write(100) ngridg
write(100) ngvec
write(100) ndmag
write(100) nspinor
write(100) fsmtype
write(100) ftmtype
write(100) dftu
write(100) lmmaxdm
write(100) xcgrad
! muffin-tin functions are unpacked to maintain backward compatibility
allocate(rfmt(lmmaxo,nrmtmax,natmtot))
if (spinpol) then
  allocate(rvfmt(lmmaxo,nrmtmax,natmtot,ndmag))
  allocate(rvfcmt(lmmaxo,nrcmtmax,natmtot,ndmag))
end if
! write the density
do ias=1,natmtot
  is=idxis(ias)
  call rfmtpack(.false.,nrmt(is),nrmti(is),rhomt(:,ias),rfmt(:,:,ias))
end do
write(100) rfmt,rhoir
! write the Coulomb potential
do ias=1,natmtot
  is=idxis(ias)
  call rfmtpack(.false.,nrmt(is),nrmti(is),vclmt(:,ias),rfmt(:,:,ias))
end do
write(100) rfmt,vclir
! write the exchange-correlation potential
do ias=1,natmtot
  is=idxis(ias)
  call rfmtpack(.false.,nrmt(is),nrmti(is),vxcmt(:,ias),rfmt(:,:,ias))
end do
write(100) rfmt,vxcir
! write the Kohn-Sham effective potential
do ias=1,natmtot
  is=idxis(ias)
  call rfmtpack(.false.,nrmt(is),nrmti(is),vsmt(:,ias),rfmt(:,:,ias))
end do
write(100) rfmt,vsir
if (spinpol) then
! write the magnetisation, exchange-correlation and effective magnetic fields
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      call rfmtpack(.false.,nrmt(is),nrmti(is),magmt(:,ias,idm), &
       rvfmt(:,:,ias,idm))
    end do
  end do
  write(100) rvfmt,magir
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      call rfmtpack(.false.,nrmt(is),nrmti(is),bxcmt(:,ias,idm), &
       rvfmt(:,:,ias,idm))
    end do
  end do
  write(100) rvfmt,bxcir
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      call rfmtpack(.false.,nrcmt(is),nrcmti(is),bsmt(:,ias,idm), &
       rvfcmt(:,:,ias,idm))
    end do
  end do
  write(100) rvfcmt,bsir
! write fixed spin moment magnetic fields
  if (fsmtype.ne.0) then
    write(100) bfsmc
    write(100) bfsmcmt
  end if
end if
! write the meta-GGA exchange-correlation potential
if (xcgrad.eq.4) then
  do ias=1,natmtot
    is=idxis(ias)
    call rfmtpack(.false.,nrmt(is),nrmti(is),wxcmt(:,ias),rfmt(:,:,ias))
  end do
  write(100) rfmt,wxcir
end if
! write the potential matrix in each muffin-tin
if ((dftu.ne.0).or.(ftmtype.ne.0)) then
  write(100) vmatmt
end if
! write the fixed tensor moment potential matrix
if (ftmtype.ne.0) then
  write(100) vmftm
end if
close(100)
deallocate(rfmt)
if (spinpol) deallocate(rvfmt,rvfcmt)
end subroutine
!EOC

