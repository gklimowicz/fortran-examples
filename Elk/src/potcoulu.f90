
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potcoulu
use modmain
use modulr
use modmpi
use modomp
implicit none
! local variables
integer ifq,is,ias
integer nr,nri,ir
integer nrc,nrci,i
integer n,lp,nthd
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:),zfmt(:)
call holdthd(nfqrz/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zvclmt,zfmt) &
!$OMP PRIVATE(ias,is,nr,nri) &
!$OMP PRIVATE(i,ir,nrc,nrci) &
!$OMP NUM_THREADS(nthd)
allocate(zrhomt(npmtmax,natmtot),zrhoir(ngtot))
allocate(zvclmt(npmtmax,natmtot),zfmt(npcmtmax))
!$OMP DO
do ifq=1,nfqrz
! distribute among MPI processes
  if (mod(ifq-1,np_mpi).ne.lp_mpi) cycle
! convert the complex muffin-tin density from coarse to fine radial mesh
  do ias=1,natmtot
    is=idxis(ias)
    zrhomt(1:npcmt(is),ias)=rhoqmt(1:npcmt(is),ias,ifq)
  end do
  call zfmtctof(zrhomt)
! solve the complex Poisson's equation in the muffin-tins
  call genzvclmt(nrmt,nrmti,nrmtmax,rlmt,wprmt,npmtmax,zrhomt,zvclmt)
! add the nuclear monopole potentials for Q=0
  if (ifq.eq.1) then
    do ias=1,natmtot
      is=idxis(ias)
      nr=nrmt(is)
      nri=nrmti(is)
      i=1
      do ir=1,nri
        zvclmt(i,ias)=zvclmt(i,ias)+vcln(ir,is)
        i=i+lmmaxi
      end do
      do ir=nri+1,nr
        zvclmt(i,ias)=zvclmt(i,ias)+vcln(ir,is)
        i=i+lmmaxo
      end do
    end do
  end if
! convert the interstitial density from coarse to fine grid
  call zfirctof(rhoqir(:,ifq),zrhoir)
! solve Poisson's equation in the entire unit cell
  call zpotcoul(nrmt,nrmti,npmt,nrmtmax,rlmt,ngridg,igfft,ngvec,gqc(:,ifq), &
   gclgq(:,ifq),ngvec,jlgqrmt(:,:,:,ifq),ylmgq(:,:,ifq),sfacgq(:,:,ifq),zrhoir,&
   npmtmax,zvclmt,vsqir(:,ifq))
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
! convert from fine to coarse radial mesh
    call zfmtftoc(nrc,nrci,zvclmt(:,ias),zfmt)
! convert to spherical coordinates
    call zbsht(nrc,nrci,zfmt,vsqmt(:,ias,ifq))
  end do
end do
!$OMP END DO
deallocate(zrhomt,zrhoir,zvclmt,zfmt)
!$OMP END PARALLEL
call freethd(nthd)
! broadcast potentials to every MPI process
if (np_mpi.gt.1) then
  n=npcmtmax*natmtot
  do ifq=1,nfqrz
    lp=mod(ifq-1,np_mpi)
    call mpi_bcast(vsqmt(:,:,ifq),n,mpi_double_complex,lp,mpicom,ierror)
  end do
  do ifq=1,nfqrz
    lp=mod(ifq-1,np_mpi)
    call mpi_bcast(vsqir(:,ifq),ngtot,mpi_double_complex,lp,mpicom,ierror)
  end do
end if
end subroutine

