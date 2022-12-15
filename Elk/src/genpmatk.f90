
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpmatk
! !INTERFACE:
subroutine genpmatk(ngp,igpig,vgpc,wfmt,wfgp,pmat)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp   : number of G+p-vectors (in,integer(nspnfv))
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
!   vgpc  : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax,nspnfv))
!   wfmt  : muffin-tin wavefunction in spherical harmonics
!           (in,complex(npcmtmax,natmtot,nspinor,nstsv))
!   wfgp  : interstitial wavefunction in plane wave basis
!           (in,complex(ngkmax,nspinor,nstsv))
!   pmat  : momentum matrix elements (out,complex(nstsv,nstsv,3))
! !DESCRIPTION:
!   Calculates the momentum matrix elements
!   $$ P_{ij}=\int d^3r\,\Psi_{i{\bf k}}^*({\bf r})\left(-i\nabla
!    +\frac{1}{4c^2}\left[\vec{\sigma}\times\nabla V_s({\bf r})\right]\right)
!    \Psi_{j{\bf k}}({\bf r}), $$
!   where $V_s$ is the Kohn-Sham effective potential. The second term in the
!   brackets is only calculated if spin-orbit coupling is enabled. See Rathgen
!   and Katsnelson, {\it Physica Scripta} {\bf T109}, 170 (2004).
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!   Fixed bug found by Juergen Spitaler, September 2006 (JKD)
!   Added spin-orbit correction, July 2010 (JKD)
!   Fixed bug found by Koichi Kitahara, January 2014 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
real(8), intent(in) :: vgpc(3,ngkmax,nspnfv)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfgp(ngkmax,nspinor,nstsv)
complex(8), intent(out) :: pmat(nstsv,nstsv,3)
! local variables
integer ist,jst,ispn,jspn
integer is,ia,ias
integer nrc,nrci,npc
integer igp,ifg,i
real(8) cso
complex(8) z1,z2,z11,z12,z21,z22,z31,z32
! automatic arrays
real(8) rfmt(npcmtmax)
complex(8) gwfmt(npcmtmax,3,nspinor),gvmt(npcmtmax,3)
complex(8) zfmt1(npcmtmax,nspinor),zfmt2(npcmtmax,3,nspinor)
complex(8) z(ngkmax)
! allocatable arrays
complex(8), allocatable :: gwfir(:,:)
! external functions
complex(8), external :: zdotc
! coefficient of spin-orbit coupling
cso=1.d0/(4.d0*solsc**2)
! zero the momentum matrix elements array
pmat(:,:,:)=0.d0
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute gradient of potential for spin-orbit correction if required
    if (spinorb) then
      call rfmtftoc(nrc,nrci,vsmt(:,ias),rfmt)
      call rtozfmt(nrc,nrci,rfmt,zfmt1)
      call gradzfmt(nrc,nrci,rlcmt(:,-1,is),wcrcmt(:,:,is),zfmt1,npcmtmax,gvmt)
! convert to spherical coordinates
      do i=1,3
        zfmt1(1:npc,1)=gvmt(1:npc,i)
        call zbsht(nrc,nrci,zfmt1,gvmt(:,i))
      end do
    end if
    do jst=1,nstsv
      do ispn=1,nspinor
! compute the gradient of the wavefunction
        call gradzfmt(nrc,nrci,rlcmt(:,-1,is),wcrcmt(:,:,is), &
         wfmt(:,ias,ispn,jst),npcmtmax,gwfmt(:,:,ispn))
      end do
! add spin-orbit correction if required
      if (spinorb) then
        do ispn=1,nspinor
! convert wavefunction to spherical coordinates
          call zbsht(nrc,nrci,wfmt(:,ias,ispn,jst),zfmt1(:,ispn))
        end do
! compute i σ x (grad V(r)) φ(r)
        do i=1,npc
          z1=zfmt1(i,1)
          z1=cmplx(-aimag(z1),dble(z1),8)
          z2=zfmt1(i,2)
          z2=cmplx(-aimag(z2),dble(z2),8)
          z11=gvmt(i,1)*z1; z12=gvmt(i,1)*z2
          z21=gvmt(i,2)*z1; z22=gvmt(i,2)*z2
          z31=gvmt(i,3)*z1; z32=gvmt(i,3)*z2
          zfmt2(i,1,1)=cmplx(aimag(z32),-dble(z32),8)-z21
          zfmt2(i,1,2)=cmplx(-aimag(z31),dble(z31),8)+z22
          zfmt2(i,2,1)=z11-z32
          zfmt2(i,2,2)=-z12-z31
          zfmt2(i,3,1)=cmplx(-aimag(z12),dble(z12),8)+z22
          zfmt2(i,3,2)=cmplx(aimag(z11),-dble(z11),8)+z21
        end do
! convert to spherical harmonics and add to wavefunction gradient
        do ispn=1,nspinor
          do i=1,3
            call zfsht(nrc,nrci,zfmt2(:,i,ispn),zfmt1)
            gwfmt(1:npc,i,ispn)=gwfmt(1:npc,i,ispn)+cso*zfmt1(1:npc,1)
          end do
        end do
      end if
      do i=1,3
        do ispn=1,nspinor
! apply the radial integral weights
          call zfmtwr(nrc,nrci,wrcmt(:,is),gwfmt(:,i,ispn))
! compute the overlaps
          do ist=1,jst
            pmat(ist,jst,i)=pmat(ist,jst,i)+zdotc(npc,wfmt(:,ias,ispn,ist),1, &
             gwfmt(:,i,ispn),1)
          end do
        end do
      end do
    end do
  end do
end do
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
allocate(gwfir(ngtc,3))
do jst=1,nstsv
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! compute the gradient
    gwfir(:,:)=0.d0
    do igp=1,ngp(jspn)
      ifg=igfc(igpig(igp,jspn))
      z1=wfgp(igp,ispn,jst)
      gwfir(ifg,1:3)=vgpc(1:3,igp,jspn)*cmplx(-aimag(z1),dble(z1),8)
    end do
    do i=1,3
! Fourier transform to real-space
      call zfftifc(3,ngdgc,1,gwfir(:,i))
! multiply by coarse characteristic function
      gwfir(:,i)=gwfir(:,i)*cfrc(:)
! Fourier transform back to G-space
      call zfftifc(3,ngdgc,-1,gwfir(:,i))
    end do
! find the overlaps
    do i=1,3
      do igp=1,ngp(jspn)
        ifg=igfc(igpig(igp,jspn))
        z(igp)=gwfir(ifg,i)
      end do
      do ist=1,jst
        pmat(ist,jst,i)=pmat(ist,jst,i)+zdotc(ngp(jspn),wfgp(:,ispn,ist),1,z,1)
      end do
    end do
  end do
end do
deallocate(gwfir)
! multiply by -i and set lower triangular part
do i=1,3
  do jst=1,nstsv
    do ist=1,jst-1
      z1=pmat(ist,jst,i)
      z1=cmplx(aimag(z1),-dble(z1),8)
      pmat(ist,jst,i)=z1
      pmat(jst,ist,i)=conjg(z1)
    end do
    pmat(jst,jst,i)=aimag(pmat(jst,jst,i))
  end do
end do
end subroutine
!EOC

