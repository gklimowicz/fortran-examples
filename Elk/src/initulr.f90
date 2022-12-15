
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initulr
use modmain
use modulr
use modomp
implicit none
! local variables
integer ik0,ik,ist,jst
integer iq,jq,ifq,ig
integer n,i,i1,i2,i3,nthd
real(8) t1
! allocatable arrays
integer, allocatable :: idx(:)
real(8), allocatable :: jlgqr(:,:)
! allocate long-range density and magnetisation arrays
if (allocated(rhormt)) deallocate(rhormt)
allocate(rhormt(npcmtmax,natmtot,nqpt))
if (allocated(rhorir)) deallocate(rhorir)
allocate(rhorir(ngtc,nqpt))
if (allocated(magrmt)) deallocate(magrmt)
if (allocated(magrir)) deallocate(magrir)
if (spinpol) then
  allocate(magrmt(npcmtmax,natmtot,ndmag,nqpt))
  allocate(magrir(ngtc,ndmag,nqpt))
end if
if (allocated(rhoqmt)) deallocate(rhoqmt)
allocate(rhoqmt(npcmtmax,natmtot,nfqrz))
if (allocated(rhoqir)) deallocate(rhoqir)
allocate(rhoqir(ngtc,nfqrz))
if (allocated(chgmtru)) deallocate(chgmtru)
allocate(chgmtru(natmtot,nqpt))
if (allocated(magqmt)) deallocate(magqmt)
if (allocated(magqir)) deallocate(magqir)
if (allocated(mommtru)) deallocate(mommtru)
if (spinpol) then
  allocate(magqmt(npcmtmax,natmtot,ndmag,nfqrz))
  allocate(magqir(ngtc,ndmag,nfqrz))
  allocate(mommtru(ndmag,natmtot,nqpt))
end if
! allocate potential and magnetic field arrays
if (allocated(vclq)) deallocate(vclq)
allocate(vclq(nfqrz))
if (allocated(bfcq)) deallocate(bfcq)
if (allocated(bfcmtq)) deallocate(bfcmtq)
if (spinpol) then
  allocate(bfcq(ndmag,nfqrz))
  allocate(bfcmtq(natmtot,ndmag,nfqrz))
end if
! combined target array for Kohn-Sham potential and magnetic field
if (allocated(vsbsq)) deallocate(vsbsq)
n=(npcmtmax*natmtot+ngtot)*nfqrz
if (spinpol) n=n*(1+ndmag)
allocate(vsbsq(n))
! zero the array
vsbsq(:)=0.d0
! associate pointer arrays with target
vsqmt(1:npcmtmax,1:natmtot,1:nfqrz)=>vsbsq(1:)
i=npcmtmax*natmtot*nfqrz+1
vsqir(1:ngtot,1:nfqrz)=>vsbsq(i:)
if (spinpol) then
  i=i+ngtot*nfqrz
  bsqmt(1:npcmtmax,1:natmtot,1:ndmag,1:nfqrz)=>vsbsq(i:)
  i=i+npcmtmax*natmtot*ndmag*nfqrz
  bsqir(1:ngtot,1:ndmag,1:nfqrz)=>vsbsq(i:)
end if
! G+Q-vector arrays
if (allocated(vgqc)) deallocate(vgqc)
allocate(vgqc(3,ngvec,nfqrz))
if (allocated(gqc)) deallocate(gqc)
allocate(gqc(ngvec,nfqrz))
if (allocated(ylmgq)) deallocate(ylmgq)
allocate(ylmgq(lmmaxo,ngvec,nfqrz))
if (allocated(sfacgq)) deallocate(sfacgq)
allocate(sfacgq(ngvec,natmtot,nfqrz))
if (allocated(gclq)) deallocate(gclq)
allocate(gclq(nqpt))
if (allocated(gclgq)) deallocate(gclgq)
allocate(gclgq(ngvec,nfqrz))
if (allocated(jlgqrmt)) deallocate(jlgqrmt)
allocate(jlgqrmt(0:lnpsd,ngvec,nspecies,nfqrz))
if (allocated(expqmt)) deallocate(expqmt)
allocate(expqmt(npcmtmax,natmtot,nqpt))
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(jlgqr,iq,ig,t1) &
!$OMP PRIVATE(i1,i2,i3,jq) &
!$OMP NUM_THREADS(nthd)
allocate(jlgqr(njcmax,nspecies))
!$OMP DO
do ifq=1,nfqrz
  iq=iqrzf(ifq)
  do ig=1,ngvec
! determine the G+Q-vectors
    vgqc(:,ig,ifq)=vgc(:,ig)+vqc(:,iq)
! G+Q-vector length
    gqc(ig,ifq)=sqrt(vgqc(1,ig,ifq)**2+vgqc(2,ig,ifq)**2+vgqc(3,ig,ifq)**2)
! spherical harmonics for G+Q-vectors
    call genylmv(lmaxo,vgqc(:,ig,ifq),ylmgq(:,ig,ifq))
  end do
! generate the spherical Bessel functions j_l(|G+Q|r)
  call genjlgpr(1,gqc(1,ifq),jlgqr)
! structure factors for G+Q-vectors
  call gensfacgp(ngvec,vgqc(:,:,ifq),ngvec,sfacgq(:,:,ifq))
! generate the Coulomb Green's function in Q-space with small Q cut-off
  t1=sqrt(vqc(1,iq)**2+vqc(2,iq)**2+vqc(3,iq)**2)
  if (t1.gt.q0cut+epslat) then
    gclq(iq)=fourpi/t1**2
  else
    gclq(iq)=0.d0
  end if
! generate the Coulomb Green's function in G+Q-space
  call gengclgq(.true.,iq,ngvec,gqc(:,ifq),gclgq(:,ifq))
! compute the spherical Bessel functions j_l(|G+Q|R_mt)
  call genjlgprmt(lnpsd,ngvec,gqc(:,ifq),ngvec,jlgqrmt(:,:,:,ifq))
! generate phase factor functions exp(iQ.r) in each muffin-tin
  call genexpmt(1,jlgqr,ylmgq(:,:,ifq),ngvec,sfacgq(:,:,ifq),expqmt(:,:,iq))
! store the phase factor function for -Q
  i1=-ivq(1,iq); i2=-ivq(2,iq); i3=-ivq(3,iq)
  if ((i1.ge.intq(1,1)).and.(i1.le.intq(2,1)).and. &
      (i2.ge.intq(1,2)).and.(i2.le.intq(2,2)).and. &
      (i3.ge.intq(1,3)).and.(i3.le.intq(2,3)).and.(ifq.gt.1)) then
    jq=ivqiq(i1,i2,i3)
    expqmt(:,:,jq)=conjg(expqmt(:,:,iq))
  end if
end do
!$OMP END DO
deallocate(jlgqr)
!$OMP END PARALLEL
call freethd(nthd)
! number of long-range states
nstulr=nstsv*nkpa
! allocate eigenvalue array
if (allocated(evalu)) deallocate(evalu)
allocate(evalu(nstulr,nkpt0))
! allocate the occupation number array
if (allocated(occulr)) deallocate(occulr)
allocate(occulr(nstulr,nkpt0))
! initialise the occupation numbers
allocate(idx(nstulr))
do ik0=1,nkpt0
  ik=(ik0-1)*nkpa+1
  call sortidx(nstulr,occsv(1,ik),idx)
  do ist=1,nstulr
    i=idx(nstulr-ist+1)-1
    ik=(ik0-1)*nkpa+i/nstsv+1
    jst=mod(i,nstsv)+1
    occulr(ist,ik0)=occsv(jst,ik)
  end do
end do
deallocate(idx)
! zero the timing variables
timemat=0.d0
timesv=0.d0
timerho=0.d0
timepot=0.d0
end subroutine

