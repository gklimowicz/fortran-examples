
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vmatmtdu
use modmain
use moddftu
use modtest
implicit none
! local variables
integer ispn,jspn
integer is,ia,ias,idu
integer l,m1,m2,m3,m4
integer ll,nm,lma,lmb
integer lm1,lm2,lm3,lm4
real(8) u,j,n,n0
real(8) mg(3),mg0(3)
real(8) v,edc,sm
complex(8) zsm,z1,z2
! automatic arrays
real(8) vee(-lmaxdm:lmaxdm,-lmaxdm:lmaxdm,-lmaxdm:lmaxdm,-lmaxdm:lmaxdm)
complex(8) dm(lmmaxdm,nspinor,lmmaxdm,nspinor)
complex(8) dms(nspinor,nspinor)
! zero the DFT+U energy for each atom
engyadu(:,:)=0.d0
do idu=1,ndftu
  is=isldu(1,idu)
  l=isldu(2,idu)
  ll=l*(l+1)+1
  nm=2*l+1
  lma=l**2+1; lmb=lma+2*l
! calculate the Coulomb matrix elements
  call genveedu(idu,u,j,vee)
  if ((abs(u).lt.1.d-10).and.(abs(j).lt.1.d-10)) cycle
! begin loop over atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! copy the density matrix for this atom
    dm(:,:,:,:)=dmatmt(:,:,:,:,ias)
! spin-unpolarised: scale density matrix so that it represents one spin channel
! (thanks to Mike Bruckhoff for this)
    if (.not.spinpol) dm(:,:,:,:)=0.5d0*dm(:,:,:,:)
! trace of density matrix for each spin
    dms(:,:)=0.d0
    do ispn=1,nspinor
      do jspn=1,nspinor
        zsm=0.d0
        do lm1=lma,lmb
          zsm=zsm+dm(lm1,ispn,lm1,jspn)
        end do
        dms(ispn,jspn)=dms(ispn,jspn)+zsm
      end do
    end do
! trace over spin
    n=dble(dms(1,1))
    if (spinpol) n=n+dble(dms(2,2))
    n0=n/dble(nspinor*nm)
! magnetisation
    if (spinpol) then
      mg(:)=0.d0
      mg(3)=dble(dms(1,1)-dms(2,2))
! non-collinear terms
      if (ncmag) then
        mg(1)=dble(dms(1,2)+dms(2,1))
        mg(2)=dble(zi*(dms(1,2)-dms(2,1)))
      end if
      mg0(:)=mg(:)/dble(nspinor*nm)
    end if
!---------------------------------!
!     around mean field (AFM)     !
!---------------------------------!
    if (dftu.eq.2) then
! modify density matrices
      do lm1=lma,lmb
        if (spinpol) then
          dm(lm1,1,lm1,1)=dm(lm1,1,lm1,1)-(n0+mg0(3))
          dm(lm1,2,lm1,2)=dm(lm1,2,lm1,2)-(n0-mg0(3))
! non-collinear terms
          if (ncmag) then
            dm(lm1,1,lm1,2)=dm(lm1,1,lm1,2)-(mg0(1)-zi*mg0(2))
            dm(lm1,2,lm1,1)=dm(lm1,2,lm1,1)-(mg0(1)+zi*mg0(2))
          end if
        else
! spin-unpolarised case
          dm(lm1,1,lm1,1)=dm(lm1,1,lm1,1)-n0
        end if
      end do
    end if
!------------------------------------!
!     DFT+U potential and energy     !
!------------------------------------!
! begin loops over m1 and m2
    do m1=-l,l
      lm1=ll+m1
      do m2=-l,l
        lm2=ll+m2
! begin loops over m3 and m4
        do m3=-l,l
          lm3=ll+m3
          do m4=-l,l
            lm4=ll+m4
            v=vee(m1,m3,m2,m4)
            do ispn=1,nspinor
              do jspn=1,nspinor
                z1=dm(lm2,ispn,lm1,ispn)*dm(lm4,jspn,lm3,jspn)
                z2=dm(lm4,jspn,lm1,ispn)*dm(lm2,ispn,lm3,jspn)
                engyadu(ia,idu)=engyadu(ia,idu)+dble(z1-z2)*v
                vmatmt(lm1,ispn,lm2,ispn,ias)=vmatmt(lm1,ispn,lm2,ispn,ias) &
                 +dm(lm4,jspn,lm3,jspn)*v
                vmatmt(lm1,ispn,lm4,jspn,ias)=vmatmt(lm1,ispn,lm4,jspn,ias) &
                 -dm(lm2,ispn,lm3,jspn)*v
              end do
            end do
! end loops over m3 and m4
          end do
        end do
! end loops over m1 and m2
      end do
    end do
! multiply energy by factor 1/2
    engyadu(ia,idu)=0.5d0*engyadu(ia,idu)
!-----------------------------------------------------------------!
!     fully localised limit (FLL) double counting corrections     !
!-----------------------------------------------------------------!
    if (dftu.eq.1) then
      if (spinpol) then
! spin-polarised
        if (ncmag) then
! non-collinear case
! correction to the energy
          edc=0.5d0*u*n*(n-1.d0)
          edc=edc-0.5d0*j*dble(dms(1,1)*(dms(1,1)-1.d0))
          edc=edc-0.5d0*j*dble(dms(2,2)*(dms(2,2)-1.d0))
          edc=edc-0.5d0*j*dble(dms(1,2)*dms(2,1))
          edc=edc-0.5d0*j*dble(dms(2,1)*dms(1,2))
! correction to the potential
          do lm1=lma,lmb
            vmatmt(lm1,1,lm1,1,ias)=vmatmt(lm1,1,lm1,1,ias) &
             -u*(n-0.5d0)+j*(dms(1,1)-0.5d0)
            vmatmt(lm1,2,lm1,2,ias)=vmatmt(lm1,2,lm1,2,ias) &
             -u*(n-0.5d0)+j*(dms(2,2)-0.5d0)
            vmatmt(lm1,1,lm1,2,ias)=vmatmt(lm1,1,lm1,2,ias)+j*dms(1,2)
            vmatmt(lm1,2,lm1,1,ias)=vmatmt(lm1,2,lm1,1,ias)+j*dms(2,1)
          end do
        else
! collinear case
! correction to the energy
          edc=0.5d0*u*n*(n-1.d0)
          edc=edc-0.5d0*j*dble(dms(1,1)*(dms(1,1)-1.d0))
          edc=edc-0.5d0*j*dble(dms(2,2)*(dms(2,2)-1.d0))
! correction to the potential
          do lm1=lma,lmb
            vmatmt(lm1,1,lm1,1,ias)=vmatmt(lm1,1,lm1,1,ias) &
             -u*(n-0.5d0)+j*(dms(1,1)-0.5d0)
            vmatmt(lm1,2,lm1,2,ias)=vmatmt(lm1,2,lm1,2,ias) &
             -u*(n-0.5d0)+j*(dms(2,2)-0.5d0)
          end do
        end if
      else
! spin-unpolarised
! correction to the energy
        edc=0.5d0*u*n*(n-1.d0)
        edc=edc-0.5d0*j*n*(n-1.d0)
! correction to the potential
        do lm1=lma,lmb
          vmatmt(lm1,1,lm1,1,ias)=vmatmt(lm1,1,lm1,1,ias)-u*(n-0.5d0) &
           +j*(n-0.5d0)
        end do
      end if
      engyadu(ia,idu)=engyadu(ia,idu)-edc
    end if
!---------------------------------------------------------!
!     subtract DFT+U potential contribution to energy     !
!---------------------------------------------------------!
! trace of dmatmt times vmatmt
    sm=0.d0
    do ispn=1,nspinor
      do lm1=lma,lmb
        do jspn=1,nspinor
          do lm2=lma,lmb
            sm=sm+dble(dm(lm1,ispn,lm2,jspn)*vmatmt(lm2,jspn,lm1,ispn,ias))
          end do
        end do
      end do
    end do
    engyadu(ia,idu)=engyadu(ia,idu)-sm
! end loop over atoms
  end do
! end loop over species
end do
! write DFT+U energy for each atom to test file
call writetest(800,'DFT+U energy for each atom',nv=natmmax*ndftu,tol=1.d-4, &
 rva=engyadu)
! write U and J parameters to test file
call writetest(810,'U and J parameters',nv=2*ndftu,tol=1.d-4,rva=ujdu)
end subroutine

