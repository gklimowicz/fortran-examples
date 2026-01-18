
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genjprk(ik,lock)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
integer(8), intent(in) :: lock(natmtot)
! local variables
integer ispn,jspn,nst,ist,jst
integer is,ia,ias,nrc,nrci,npc
integer igk,ifg,i
real(8) wo
complex(8) z1
! automatic arrays
integer idx(nstsv)
real(8) rfmt(npcmtmax)
complex(8) gwfmt(npcmtmax,3),zfmt1(npcmtmax),zfmt2(npcmtmax)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgk(:,:,:)
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
! find the matching coefficients
do ispn=1,nspnfv
  call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
! get the eigenvectors from file
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! count and index the occupied states
nst=0
do ist=1,nstsv
  if (abs(occsv(ist,ik)).lt.epsocc) cycle
  nst=nst+1
  idx(nst)=ist
end do
! calculate the second-variational wavefunctions for occupied states
allocate(wfmt(npcmtmax,natmtot,nspinor,nst),wfgk(ngkmax,nspinor,nst))
call genwfsv(.true.,.true.,nst,idx,ngdgc,igfc,ngk(:,ik),igkig(:,:,ik),apwalm, &
 evecfv,evecsv,wfmt,ngkmax,wfgk)
deallocate(apwalm,evecfv,evecsv)
!-------------------------------------------------!
!     muffin-tin paramagnetic current density     !
!-------------------------------------------------!
do ist=1,nst
  jst=idx(ist)
  wo=wkpt(ik)*occsv(jst,ik)
  do ispn=1,nspinor
    do is=1,nspecies
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! compute the gradient of the wavefunction
        call gradzfmt(nrc,nrci,rlcmt(:,-1,is),wcrcmt(:,:,is), &
         wfmt(:,ias,ispn,ist),npcmtmax,gwfmt)
! convert wavefunction to spherical coordinates and conjugate
        call zbsht(nrc,nrci,wfmt(:,ias,ispn,ist),zfmt1)
        zfmt1(1:npc)=conjg(zfmt1(1:npc))
        do i=1,3
! convert wavefunction gradient to spherical coordinates
          call zbsht(nrc,nrci,gwfmt(:,i),zfmt2)
! compute the partial current density
          rfmt(1:npc)=aimag(zfmt1(1:npc)*zfmt2(1:npc))
          call omp_set_lock(lock(ias))
          jrmt(1:npc,ias,i)=jrmt(1:npc,ias,i)+wo*rfmt(1:npc)
          call omp_unset_lock(lock(ias))
        end do
      end do
    end do
  end do
end do
deallocate(wfmt)
!---------------------------------------------------!
!     interstitial paramagnetic current density     !
!---------------------------------------------------!
allocate(zfft1(ngtc),zfft2(ngtc))
do ist=1,nst
  jst=idx(ist)
  wo=wkpt(ik)*occsv(jst,ik)/omega
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! Fourier transform to real-space and conjugate
    zfft1(:)=0.d0
    do igk=1,ngk(jspn,ik)
      ifg=igfc(igkig(igk,jspn,ik))
      zfft1(ifg)=wfgk(igk,ispn,ist)
    end do
    call zfftifc(3,ngdgc,1,zfft1)
    zfft1(:)=conjg(zfft1(:))
    do i=1,3
! compute the gradient of the wavefunction
      zfft2(:)=0.d0
      do igk=1,ngk(jspn,ik)
        ifg=igfc(igkig(igk,jspn,ik))
        z1=wfgk(igk,ispn,ist)
        zfft2(ifg)=vgkc(i,igk,jspn,ik)*cmplx(-aimag(z1),dble(z1),8)
      end do
      call zfftifc(3,ngdgc,1,zfft2)
!$OMP CRITICAL(currentk_)
      jrir(1:ngtc,i)=jrir(1:ngtc,i)+wo*aimag(zfft1(1:ngtc)*zfft2(1:ngtc))
!$OMP END CRITICAL(currentk_)
    end do
  end do
end do
deallocate(wfgk,zfft1,zfft2)
end subroutine

