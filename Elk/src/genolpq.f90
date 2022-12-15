
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genolpq(nst,expqmt,ngpq,igpqig,wfmt,wfir,wfmtq,wfgpq,oq)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: nst
complex(8), intent(in) :: expqmt(npcmtmax,natmtot)
integer, intent(in) :: ngpq(nspnfv),igpqig(ngkmax,nspnfv)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nst)
complex(8), intent(in) :: wfir(ngtot,nspinor,nst)
complex(8), intent(in) :: wfmtq(npcmtmax,natmtot,nspinor,nst)
complex(8), intent(in) :: wfgpq(ngkmax,nspinor,nst)
complex(8), intent(out) :: oq(nst,nst)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,nrc,nrci
integer npc,igpq,nthd
real(8) t1
! automatic arrays
complex(8) wfmt1(npcmtmax),z(ngkmax)
! allocatable arrays
complex(8), allocatable :: wfir1(:)
! external functions
complex(8), external :: zdotc
call holdthd(nst,nthd)
! zero the matrix elements
oq(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,ispn,ias,is) &
!$OMP PRIVATE(nrc,nrci,npc,ist) &
!$OMP NUM_THREADS(nthd)
do jst=1,nst
  do ispn=1,nspinor
    do ias=1,natmtot
      is=idxis(ias)
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
! multiply by local phase factor function exp(iq.r)
      wfmt1(1:npc)=expqmt(1:npc,ias)*wfmt(1:npc,ias,ispn,jst)
! apply the radial integral weights
      call zfcmtwr(nrc,nrci,wrcmt(:,is),wfmt1)
! compute the inner products
      do ist=1,nst
        oq(ist,jst)=oq(ist,jst)+zdotc(npc,wfmtq(:,ias,ispn,ist),1,wfmt1,1)
      end do
    end do
  end do
end do
!$OMP END PARALLEL DO
!---------------------------!
!     interstitial part     !
!---------------------------!
t1=sqrt(omega)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,z,ispn,jspn) &
!$OMP PRIVATE(igpq,ist) &
!$OMP NUM_THREADS(nthd)
allocate(wfir1(ngtot))
!$OMP DO
do jst=1,nst
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! multiply wavefunction by characteristic function
    wfir1(:)=wfir(:,ispn,jst)*cfunir(:)
! Fourier transform to G+p+q-space
    call zfftifc(3,ngridg,-1,wfir1)
    do igpq=1,ngpq(jspn)
      z(igpq)=wfir1(igfft(igpqig(igpq,jspn)))
    end do
    do ist=1,nst
      oq(ist,jst)=oq(ist,jst)+t1*zdotc(ngpq(jspn),wfgpq(:,ispn,ist),1,z,1)
    end do
  end do
end do
!$OMP END DO
deallocate(wfir1)
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

