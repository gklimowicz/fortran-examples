
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genhmlu(ik0,h)
use modmain
use modulr
use modomp
implicit none
! arguments
integer, intent(in) :: ik0
complex(8), intent(out) :: h(nstulr,nstulr)
! local variables
integer ik,ist,jst,ispn,nthd
integer ikpa,jkpa,iq,ifq,ngk0,igk
integer i1,i2,i3,j1,j2,j3,i,j
! automatic arrays
complex(8) vmat(nstsv,nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(4), allocatable :: wfmt(:,:,:,:),wfir(:,:,:),wfgk(:,:,:)
complex(8), allocatable :: hdb(:,:,:)
! central k-point
ik=(ik0-1)*nkpa+1
! number of G+k-vectors for central k-point
ngk0=ngk(1,ik)
! get the ground-state eigenvectors from file for central k-point
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
call getevecfv('.OUT',ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk0,vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! calculate the wavefunctions for all states of the central k-point
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfgk(ngk0,nspinor,nstsv))
call genwfsv_sp(.false.,.true.,nstsv,0,ngridg,igfft,ngk0,igkig(:,1,ik),apwalm, &
 evecfv,evecsv,wfmt,ngk0,wfgk)
deallocate(apwalm,evecfv,evecsv)
! determine the interstitial wavefunctions in real-space (without 1/sqrt(omega))
allocate(wfir(ngtot,nspinor,nstsv))
call holdthd(nstsv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ispn,igk) &
!$OMP NUM_THREADS(nthd)
do ist=1,nstsv
  do ispn=1,nspinor
    wfir(:,ispn,ist)=0.e0
    do igk=1,ngk0
      wfir(igfft(igkig(igk,1,ik)),ispn,ist)=wfgk(igk,ispn,ist)
    end do
    call cfftifc(3,ngridg,1,wfir(:,ispn,ist))
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! generate the matrix elements for all Q-vectors
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(vmat,iq,i,j,ikpa,jkpa) &
!$OMP PRIVATE(j1,j2,j3,ist,jst,i1,i2,i3) &
!$OMP NUM_THREADS(nthd)
do ifq=1,nfqrz
  iq=iqrzf(ifq)
  if (spinpol) then
    call genzvbmatk(vsqmt(:,:,ifq),vsqir(:,ifq),bsqmt(:,:,:,ifq), &
     bsqir(:,:,ifq),ngk0,igkig(:,1,ik),wfmt,wfir,wfgk,vmat)
  else
    call genzvmatk(vsqmt(:,:,ifq),vsqir(:,ifq),ngk0,igkig(:,1,ik),wfmt,wfir, &
     wfgk,vmat)
  end if
  j=0
  do jkpa=1,nkpa
    j1=ivq(1,jkpa); j2=ivq(2,jkpa); j3=ivq(3,jkpa)
    do jst=1,nstsv
      j=j+1
      do ikpa=1,jkpa-1
        i=(ikpa-1)*nstsv+1
        i1=ivq(1,ikpa)-j1; i2=ivq(2,ikpa)-j2; i3=ivq(3,ikpa)-j3
        if (ivqiq(i1,i2,i3).eq.iq) then
! copy matrix elements for kappa_i - kappa_j in Q-point set
          call zcopy(nstsv,vmat(:,jst),1,h(i,j),1)
        else if (ivqiq(-i1,-i2,-i3).eq.iq) then
! otherwise use conjugate transpose
          do ist=1,nstsv
            h(i,j)=conjg(vmat(jst,ist))
            i=i+1
          end do
        end if
      end do
! copy only the upper triangular part for Q=0
      if (ifq.eq.1) then
        i=(jkpa-1)*nstsv+1
        call zcopy(jst,vmat(:,jst),1,h(i,j),1)
      end if
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
deallocate(wfmt,wfir,wfgk)
! add the second-variational eigenvalues of k+kappa to the diagonal but in the
! basis of the states at k
do ist=1,nstsv
  h(ist,ist)=h(ist,ist)+evalsv(ist,ik)
end do
allocate(hdb(nstsv,nstsv,2:nkpa))
call gethdbulr(ik0,hdb)
do ikpa=2,nkpa
  i=(ikpa-1)*nstsv
  do jst=1,nstsv
    j=i+jst
    do ist=1,jst
      h(i+ist,j)=h(i+ist,j)+hdb(ist,jst,ikpa)
    end do
  end do
end do
deallocate(hdb)
end subroutine

