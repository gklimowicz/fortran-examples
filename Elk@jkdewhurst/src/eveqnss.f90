
! Copyright (C) 2006 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnss(ngp,igpig,apwalm,evalfv,evecfv,evalsvp,evecsv)
use modmain
use moddftu
use modomp
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
real(8), intent(in) :: evalfv(nstfv,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
real(8), intent(out) :: evalsvp(nstsv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ist,jst,ispn
integer is,ias,i,j,k
integer nrc,nrci,nrco
integer l,lm,nm,npc,npci
integer igp,ld,nthd
real(8) ts0,ts1
complex(8) zq
! automatic arrays
complex(8) wfmt2(npcmtmax,nspnfv),wfmt4(npcmtmax)
complex(8) wfmt31(npcmtmax),wfmt32(npcmtmax),wfmt33(npcmtmax)
complex(8) wfgp1(ngkmax),wfgp2(ngkmax),wfgp3(ngkmax)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:),wfir1(:,:),wfir2(:)
! external functions
complex(8), external :: zdotc
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(eveqnss): spin-unpolarised calculation")')
  write(*,*)
  stop
end if
call timesec(ts0)
ld=lmmaxdm*nspinor
call holdthd(nstfv,nthd)
! zero the second-variational Hamiltonian (stored in the eigenvector array)
evecsv(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(wfmt1(npcmtmax,nstfv,nspnfv))
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt2,wfmt31,wfmt32,wfmt33,wfmt4) &
!$OMP PRIVATE(ias,is,nrc,nrci,nrco) &
!$OMP PRIVATE(npc,npci,zq,ispn) &
!$OMP PRIVATE(ist,jst,l,nm,lm,i,j,k) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  npc=npcmt(is)
  npci=npcmti(is)
! de-phasing factor (FC, FB & LN)
  if (ssdph) zq=zqss(ias)
! compute the first-variational wavefunctions
  do ispn=1,nspnfv
    if (ssdph.and.(ispn.eq.2)) zq=conjg(zq)
!$OMP DO
    do ist=1,nstfv
      call wfmtfv(ias,ngp(ispn),apwalm(:,:,:,ias,ispn),evecfv(:,ist,ispn), &
       wfmt1(:,ist,ispn))
! de-phase if required
      if (ssdph) wfmt1(1:npc,ist,ispn)=zq*wfmt1(1:npc,ist,ispn)
    end do
!$OMP END DO
  end do
!$OMP DO
  do jst=1,nstfv
! convert wavefunction to spherical coordinates
    do ispn=1,nspnfv
      call zbsht(nrc,nrci,wfmt1(:,jst,ispn),wfmt2(:,ispn))
    end do
! apply effective magnetic field and convert to spherical harmonics
    wfmt4(1:npc)=bsmt(1:npc,ias,3)*wfmt2(1:npc,1)
    call zfsht(nrc,nrci,wfmt4,wfmt31)
    wfmt4(1:npc)=-bsmt(1:npc,ias,3)*wfmt2(1:npc,2)
    call zfsht(nrc,nrci,wfmt4,wfmt32)
    wfmt4(1:npc)=cmplx(bsmt(1:npc,ias,1),-bsmt(1:npc,ias,2),8)*wfmt2(1:npc,2)
    call zfsht(nrc,nrci,wfmt4,wfmt33)
! apply muffin-tin DFT+U potential matrix if required
    if (tvmatmt) then
      do l=0,lmaxdm
        if (tvmmt(l,ias)) then
          nm=2*l+1
          lm=l**2+1
          i=npci+lm
          if (l.le.lmaxi) then
            call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,1,lm,1,ias),ld, &
             wfmt1(lm,jst,1),lmmaxi,zone,wfmt31(lm),lmmaxi)
          end if
          call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,1,lm,1,ias),ld, &
           wfmt1(i,jst,1),lmmaxo,zone,wfmt31(i),lmmaxo)
          if (l.le.lmaxi) then
            call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,2,lm,2,ias),ld, &
             wfmt1(lm,jst,2),lmmaxi,zone,wfmt32(lm),lmmaxi)
          end if
          call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,2,lm,2,ias),ld, &
           wfmt1(i,jst,2),lmmaxo,zone,wfmt32(i),lmmaxo)
          if (l.le.lmaxi) then
            call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,1,lm,2,ias),ld, &
             wfmt1(lm,jst,2),lmmaxi,zone,wfmt33(lm),lmmaxi)
          end if
          call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,1,lm,2,ias),ld, &
           wfmt1(i,jst,2),lmmaxo,zone,wfmt33(i),lmmaxo)
        end if
      end do
    end if
! add to second-variational Hamiltonian matrix
! upper diagonal block
    call zfmtwr(nrc,nrci,wrcmt(:,is),wfmt31)
    do ist=1,jst
      evecsv(ist,jst)=evecsv(ist,jst)+zdotc(npc,wfmt1(:,ist,1),1,wfmt31,1)
    end do
    j=jst+nstfv
! lower diagonal block
    call zfmtwr(nrc,nrci,wrcmt(:,is),wfmt32)
    do ist=1,jst
      i=ist+nstfv
      evecsv(i,j)=evecsv(i,j)+zdotc(npc,wfmt1(:,ist,2),1,wfmt32,1)
    end do
! off-diagonal block
    call zfmtwr(nrc,nrci,wrcmt(:,is),wfmt33)
    do ist=1,nstfv
      evecsv(ist,j)=evecsv(ist,j)+zdotc(npc,wfmt1(:,ist,1),1,wfmt33,1)
    end do
  end do
!$OMP END DO
! end loop over atoms
end do
!$OMP END PARALLEL
deallocate(wfmt1)
!---------------------------!
!     interstitial part     !
!---------------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,wfir2,wfgp1,wfgp2,wfgp3) &
!$OMP PRIVATE(ispn,igp,ist,i,j) &
!$OMP NUM_THREADS(nthd)
allocate(wfir1(ngtot,nspnfv),wfir2(ngtot))
! begin loop over states
!$OMP DO
do jst=1,nstfv
  do ispn=1,nspnfv
    wfir1(:,ispn)=0.d0
    do igp=1,ngp(ispn)
      wfir1(igfft(igpig(igp,ispn)),ispn)=evecfv(igp,jst,ispn)
    end do
! Fourier transform wavefunction to real-space
    call zfftifc(3,ngridg,1,wfir1(:,ispn))
  end do
! multiply with magnetic field and transform to G-space
  wfir2(1:ngtot)=bsir(1:ngtot,3)*wfir1(1:ngtot,1)
  call zfftifc(3,ngridg,-1,wfir2)
  do igp=1,ngp(1)
    wfgp1(igp)=wfir2(igfft(igpig(igp,1)))
  end do
  wfir2(1:ngtot)=-bsir(1:ngtot,3)*wfir1(1:ngtot,2)
  call zfftifc(3,ngridg,-1,wfir2)
  do igp=1,ngp(2)
    wfgp2(igp)=wfir2(igfft(igpig(igp,2)))
  end do
  wfir2(1:ngtot)=cmplx(bsir(1:ngtot,1),-bsir(1:ngtot,2),8)*wfir1(1:ngtot,2)
  call zfftifc(3,ngridg,-1,wfir2)
  do igp=1,ngp(1)
    wfgp3(igp)=wfir2(igfft(igpig(igp,1)))
  end do
! add to second-variational Hamiltonian matrix
! upper diagonal block
  do ist=1,jst
    evecsv(ist,jst)=evecsv(ist,jst)+zdotc(ngp(1),evecfv(:,ist,1),1,wfgp1,1)
  end do
! lower diagonal block
  j=jst+nstfv
  do ist=1,jst
    i=ist+nstfv
    evecsv(i,j)=evecsv(i,j)+zdotc(ngp(2),evecfv(:,ist,2),1,wfgp2,1)
  end do
! off-diagonal block
  do ist=1,nstfv
    evecsv(ist,j)=evecsv(ist,j)+zdotc(ngp(1),evecfv(:,ist,1),1,wfgp3,1)
  end do
end do
!$OMP END DO
deallocate(wfir1,wfir2)
!$OMP END PARALLEL
call freethd(nthd)
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist,ispn)
  end do
end do
! diagonalise the second-variational Hamiltonian
call eveqnzh(nstsv,nstsv,evecsv,evalsvp)
call timesec(ts1)
!$OMP ATOMIC
timesv=timesv+ts1-ts0
end subroutine

