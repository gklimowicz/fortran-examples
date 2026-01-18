
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potxcu
use modmain
use modulr
use modmpi
use modomp
implicit none
! local variables
integer ifq,idm,is,ias
integer ir,npc,n,lp,nthd
complex(8) z1,z2
! allocatable arrays
real(8), allocatable :: vxcrmt(:,:,:),vxcrir(:,:)
real(8), allocatable :: bxcrmt(:,:,:,:),bxcrir(:,:,:)
real(8), allocatable :: rhomt_(:,:),rhoir_(:)
real(8), allocatable :: magmt_(:,:,:),magir_(:,:)
real(8), allocatable :: vxcmt_(:,:),bxcmt_(:,:,:)
complex(8), allocatable :: vxcqmt(:,:,:),vxcqir(:,:)
allocate(vxcrmt(npcmtmax,natmtot,nqpt),vxcrir(ngtot,nqpt))
if (spinpol) then
  allocate(bxcrmt(npcmtmax,natmtot,ndmag,nqpt))
  allocate(bxcrir(ngtot,ndmag,nqpt))
end if
call holdthd(nqpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rhomt_,rhoir_,vxcmt_) &
!$OMP PRIVATE(magmt_,magir_,bxcmt_) &
!$OMP PRIVATE(ias,is,idm) &
!$OMP NUM_THREADS(nthd)
allocate(rhomt_(npmtmax,natmtot),rhoir_(ngtot))
allocate(vxcmt_(npmtmax,natmtot))
if (spinpol) then
  allocate(magmt_(npmtmax,natmtot,ndmag),magir_(ngtot,ndmag))
  allocate(bxcmt_(npmtmax,natmtot,ndmag))
end if
!$OMP DO
do ir=1,nqpt
! distribute among MPI processes
  if (mod(ir-1,np_mpi).ne.lp_mpi) cycle
! convert muffin-tin density and magnetisation from coarse to fine radial mesh
  do ias=1,natmtot
    is=idxis(ias)
    rhomt_(1:npcmt(is),ias)=rhormt(1:npcmt(is),ias,ir)
  end do
  call rfmtctof(rhomt_)
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      magmt_(1:npcmt(is),ias,idm)=magrmt(1:npcmt(is),ias,idm,ir)
    end do
    call rfmtctof(magmt_(:,:,idm))
  end do
! convert interstitial density and magnetisation from coarse to fine grid
  call rfirctof(rhorir(:,ir),rhoir_)
  do idm=1,ndmag
    call rfirctof(magrir(:,idm,ir),magir_(:,idm))
  end do
! calculate the exchange-correlation potential and magnetic field
  call potxc(.false.,xctype,rhomt_,rhoir_,magmt_,magir_,taumt,tauir,exmt,exir, &
   ecmt,ecir,vxcmt_,vxcrir(:,ir),bxcmt_,bxcrir(:,:,ir),wxcmt,wxcir)
! convert muffin-tin potential and field from fine to coarse radial mesh
  do ias=1,natmtot
    is=idxis(ias)
    call rfmtftoc(nrcmt(is),nrcmti(is),vxcmt_(:,ias),vxcrmt(:,ias,ir))
  end do
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      call rfmtftoc(nrcmt(is),nrcmti(is),bxcmt_(:,ias,idm),bxcrmt(:,ias,idm,ir))
    end do
  end do
end do
!$OMP END DO
deallocate(rhomt_,rhoir_,vxcmt_)
if (spinpol) deallocate(magmt_,magir_,bxcmt_)
!$OMP END PARALLEL
call freethd(nthd)
! broadcast potentials and fields to every MPI process
if (np_mpi.gt.1) then
  n=npcmtmax*natmtot
  do ir=1,nqpt
    lp=mod(ir-1,np_mpi)
    call mpi_bcast(vxcrmt(:,:,ir),n,mpi_double_precision,lp,mpicom,ierror)
  end do
  do ir=1,nqpt
    lp=mod(ir-1,np_mpi)
    call mpi_bcast(vxcrir(:,ir),ngtot,mpi_double_precision,lp,mpicom,ierror)
  end do
  if (spinpol) then
    n=npcmtmax*natmtot*ndmag
    do ir=1,nqpt
      lp=mod(ir-1,np_mpi)
      call mpi_bcast(bxcrmt(:,:,:,ir),n,mpi_double_precision,lp,mpicom,ierror)
    end do
    n=ngtot*ndmag
    do ir=1,nqpt
      lp=mod(ir-1,np_mpi)
      call mpi_bcast(bxcrir(:,:,ir),n,mpi_double_precision,lp,mpicom,ierror)
    end do
  end if
end if
allocate(vxcqmt(npcmtmax,natmtot,nfqrz),vxcqir(ngtot,nfqrz))
! Fourier transform exchange-correlation potential to Q-space
call rfzfftq(-1,1,ngtot,vxcrmt,vxcrir,vxcqmt,vxcqir)
deallocate(vxcrmt,vxcrir)
! add V_xc and external Coulomb potential to Kohn-Sham potential
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(z1,ias,is,npc) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ifq=1,nfqrz
  z1=vclq(ifq)
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
    vsqmt(1:npc,ias,ifq)=vsqmt(1:npc,ias,ifq)+vxcqmt(1:npc,ias,ifq)+z1
  end do
end do
!$OMP END DO NOWAIT
!$OMP DO
do ifq=1,nfqrz
  z1=vclq(ifq)
  vsqir(:,ifq)=vsqir(:,ifq)+vxcqir(:,ifq)+z1
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
deallocate(vxcqmt,vxcqir)
if (spinpol) then
! Fourier transform the exchange-correlation magnetic field to Q-space
  call rfzfftq(-1,ndmag,ngtot,bxcrmt,bxcrir,bsqmt,bsqir)
  deallocate(bxcrmt,bxcrir)
! add external magnetic fields
  call holdthd(nfqrz,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(idm,z1,z2,ias,is,npc) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
  do ifq=1,nfqrz
    do idm=1,ndmag
      z1=bfcq(idm,ifq)
      do ias=1,natmtot
        is=idxis(ias)
        npc=npcmt(is)
        z2=z1+bfcmtq(ias,idm,ifq)
        bsqmt(1:npc,ias,idm,ifq)=bsqmt(1:npc,ias,idm,ifq)+z2
      end do
    end do
  end do
!$OMP END DO NOWAIT
!$OMP DO
  do ifq=1,nfqrz
    do idm=1,ndmag
      z1=bfcq(idm,ifq)
      bsqir(:,idm,ifq)=bsqir(:,idm,ifq)+z1
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL
  call freethd(nthd)
end if
end subroutine

