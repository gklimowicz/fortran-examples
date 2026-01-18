
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradwf2(ik,gwf2mt,gwf2ir)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(inout) :: gwf2mt(npmtmax,natmtot),gwf2ir(ngtot)
! local variables
integer ispn,jspn,nst,ist,jst
integer is,ia,ias,nrc,nrci,npc
integer igk,ifg,i
real(8) wo
complex(8) z1
! automatic arrays
integer idx(nstsv)
complex(8) gwfmt(npcmtmax,3),zfmt(npcmtmax)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgk(:,:,:)
complex(8), allocatable :: zfft(:)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
do ispn=1,nspnfv
  call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik),sfacgk(:,:,ispn,ik),&
   apwalm(:,:,:,:,ispn))
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
deallocate(apwalm)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
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
        do i=1,3
! convert gradient from spherical harmonics to spherical coordinates
          call zbsht(nrc,nrci,gwfmt(:,i),zfmt)
! add to total
          gwf2mt(1:npc,ias)=gwf2mt(1:npc,ias) &
           +wo*(dble(zfmt(1:npc))**2+aimag(zfmt(1:npc))**2)
        end do
      end do
    end do
  end do
end do
deallocate(wfmt)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(zfft(ngtc))
do ist=1,nst
  jst=idx(ist)
  wo=wkpt(ik)*occsv(jst,ik)/omega
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! compute gradient of wavefunction
    do i=1,3
      zfft(:)=0.d0
      do igk=1,ngk(jspn,ik)
        ifg=igfc(igkig(igk,jspn,ik))
        z1=wfgk(igk,ispn,ist)
        zfft(ifg)=vgkc(i,igk,jspn,ik)*cmplx(-aimag(z1),dble(z1),8)
      end do
      call zfftifc(3,ngdgc,1,zfft)
      gwf2ir(1:ngtc)=gwf2ir(1:ngtc) &
       +wo*(dble(zfft(1:ngtc))**2+aimag(zfft(1:ngtc))**2)
    end do
  end do
end do
deallocate(wfgk,zfft)
end subroutine

