
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmatk(vmt,vir,ngp,igpig,wfmt,ld,wfgp,vmat)
use modmain
use moddftu
use modomp
implicit none
! arguments
! the potential is multiplied by the radial integration weights in the
! muffin-tin and by the characteristic function in the interstitial region
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(4), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ld
complex(4), intent(in) :: wfgp(ld,nspinor,nstsv)
complex(8), intent(out) :: vmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn,n
integer is,ias,nrc,nrci,nrco
integer npc,npc2,ipco,nthd
! automatic arrays
complex(4) wfmt1(npcmtmax),wfmt2(npcmtmax),c(ngkmax)
! allocatable arrays
complex(4), allocatable :: wfir(:)
! external functions
real(4), external :: sdot
complex(4), external :: cdotc
call holdthd(nstsv,nthd)
! zero the matrix elements
vmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do ispn=1,nspinor
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    nrco=nrc-nrci
    npc=npcmt(is)
    npc2=npc*2
    ipco=npcmti(is)+1
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,wfmt2,ist) &
!$OMP NUM_THREADS(nthd)
    do jst=1,nstsv
      wfmt1(1:npc)=vmt(1:npc,ias)*wfmt(1:npc,ias,ispn,jst)
! apply muffin-tin DFT+U potential matrix if required
      if (dftu.ne.0) then
        if (any(tvmmt(0:lmaxdm,ias))) then
          call cgemm('N','N',lmmaxi,nrci,lmmaxi,cone,vmatmti(1,1,1,1,ias), &
           lmmaxi,wfmt(1,ias,ispn,jst),lmmaxi,czero,wfmt2,lmmaxi)
          call cgemm('N','N',lmmaxo,nrco,lmmaxo,cone,vmatmto(1,1,1,1,ias), &
           lmmaxo,wfmt(ipco,ias,ispn,jst),lmmaxo,czero,wfmt2(ipco),lmmaxo)
          call cfcmtwr(nrc,nrci,wrcmt(:,is),wfmt2)
          wfmt1(1:npc)=wfmt1(1:npc)+wfmt2(1:npc)
        end if
      end if
! compute the inner products
      do ist=1,jst-1
        vmat(ist,jst)=vmat(ist,jst)+cdotc(npc,wfmt(1,ias,ispn,ist),1,wfmt1,1)
      end do
      vmat(jst,jst)=vmat(jst,jst)+sdot(npc2,wfmt(1,ias,ispn,jst),1,wfmt1,1)
    end do
!$OMP END PARALLEL DO
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir,c,ispn,jspn,n,ist) &
!$OMP NUM_THREADS(nthd)
allocate(wfir(ngtot))
!$OMP DO
do jst=1,nstsv
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    n=ngp(jspn)
! Fourier transform wavefunction to real-space
    wfir(:)=0.e0
    wfir(igfft(igpig(1:n,jspn)))=wfgp(1:n,ispn,jst)
    call cfftifc(3,ngridg,1,wfir)
! apply potential to wavefunction
    wfir(1:ngtot)=vir(1:ngtot)*wfir(1:ngtot)
! Fourier transform to G+p-space
    call cfftifc(3,ngridg,-1,wfir)
    c(1:n)=wfir(igfft(igpig(1:n,jspn)))
    do ist=1,jst-1
! compute inner product
      vmat(ist,jst)=vmat(ist,jst)+cdotc(n,wfgp(1,ispn,ist),1,c,1)
    end do
    vmat(jst,jst)=vmat(jst,jst)+sdot(n*2,wfgp(1,ispn,jst),1,c,1)
  end do
end do
!$OMP END DO
deallocate(wfir)
!$OMP END PARALLEL
call freethd(nthd)
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vmat(ist,jst)=conjg(vmat(jst,ist))
  end do
end do
end subroutine

