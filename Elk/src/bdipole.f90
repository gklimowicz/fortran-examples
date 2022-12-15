
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bdipole
use modmain
implicit none
! local variables
integer idm,is,ias
integer nrc,nrci,npc
real(8) cb,t1
! automatic arrays
real(8) rfmt1(npcmtmax),rfmt2(npcmtmax)
! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:),rvfir(:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
if (.not.ncmag) then
  write(*,*)
  write(*,'("Error(bdipole): non-collinear magnetism required for inclusion of &
   &the dipole field")')
  write(*,*)
  stop
end if
! prefactor for the spin dipole magnetic field
cb=gfacte/(4.d0*solsc)
! compute the gauge invariant current density if required
if (tjr) call genjr
! allocate local arrays
allocate(rvfmt(npmtmax,natmtot,3),rvfir(ngtot,3))
allocate(zrhomt(npmtmax,natmtot),zrhoir(ngtot))
allocate(zvclmt(npmtmax,natmtot),zvclir(ngtot))
! compute the curl of the magnetisation density, i.e. the magnetisation current
call curlrvf(magmt,magir,rvfmt,rvfir)
! multiply by prefactor
rvfmt(:,:,:)=cb*rvfmt(:,:,:)
rvfir(:,:)=cb*rvfir(:,:)
! add the current density if required
if (tjr) then
  t1=1.d0/solsc
  rvfmt(:,:,:)=rvfmt(:,:,:)+t1*jrmt(:,:,:)
  rvfir(:,:)=rvfir(:,:)+t1*jrir(:,:)
end if
do idm=1,3
! transform to complex spherical harmonics
  do ias=1,natmtot
    is=idxis(ias)
    call rtozfmt(nrmt(is),nrmti(is),rvfmt(:,ias,idm),zrhomt(:,ias))
  end do
! solve Poisson's equation in the muffin-tin to find the A-field
  call genzvclmt(nrmt,nrmti,nrmtmax,rlmt,wprmt,npmtmax,zrhomt,zvclmt)
  zrhoir(:)=rvfir(:,idm)
! solve in the entire unit cell
  call zpotcoul(nrmt,nrmti,npmt,nrmtmax,rlmt,ngridg,igfft,ngvec,gc,gclg,ngvec, &
   jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)
! convert muffin-tin A-field to real spherical harmonics
  do ias=1,natmtot
    is=idxis(ias)
    call ztorfmt(nrmt(is),nrmti(is),zvclmt(:,ias),rvfmt(:,ias,idm))
  end do
! store the real part of the interstitial A-field
  rvfir(:,idm)=dble(zvclir(:))
end do
! compute the curl of A to obtain the dipole B-field
call curlrvf(rvfmt,rvfir,bdmt,bdir)
! add to the Kohn-Sham field
do idm=1,3
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
! convert to coarse radial mesh
    call rfmtftoc(nrc,nrci,bdmt(:,ias,idm),rfmt1)
! convert to spherical coordinates
    call rbsht(nrc,nrci,rfmt1,rfmt2)
    bsmt(1:npc,ias,idm)=bsmt(1:npc,ias,idm)+cb*rfmt2(1:npc)
  end do
end do
do idm=1,3
  bsir(:,idm)=bsir(:,idm)+cb*bdir(:,idm)*cfunir(:)
end do
deallocate(rvfmt,rvfir)
deallocate(zrhomt,zrhoir,zvclmt,zvclir)
end subroutine

