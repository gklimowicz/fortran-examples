
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potxc
! !INTERFACE:
subroutine potxc(tsh,xctype_,rhomt_,rhoir_,magmt_,magir_,taumt_,tauir_,exmt_, &
 exir_,ecmt_,ecir_,vxcmt_,vxcir_,bxcmt_,bxcir_,wxcmt_,wxcir_)
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Computes the exchange-correlation potential and energy density. In the
!   muffin-tin, the density is transformed from spherical harmonic coefficients
!   $\rho_{lm}$ to spherical coordinates $(\theta,\phi)$ with a backward
!   spherical harmonic transformation (SHT). Once calculated, the
!   exchange-correlation potential and energy density are transformed with a
!   forward SHT.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: xctype_(3)
real(8), intent(in) :: rhomt_(npmtmax,natmtot),rhoir_(ngtot)
real(8), intent(in) :: magmt_(npmtmax,natmtot,ndmag),magir_(ngtot,ndmag)
real(8), intent(in) :: taumt_(npmtmax,natmtot,nspinor),tauir_(ngtot,nspinor)
real(8), intent(out) :: exmt_(npmtmax,natmtot),exir_(ngtot)
real(8), intent(out) :: ecmt_(npmtmax,natmtot),ecir_(ngtot)
real(8), intent(out) :: vxcmt_(npmtmax,natmtot),vxcir_(ngtot)
real(8), intent(out) :: bxcmt_(npmtmax,natmtot,ndmag),bxcir_(ngtot,ndmag)
real(8), intent(out) :: wxcmt_(npmtmax,natmtot),wxcir_(ngtot)
! local variables
integer ias,nthd1,nthd2
call holdthd(2,nthd1)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(nthd2,ias) &
!$OMP NUM_THREADS(nthd1)
!$OMP SECTIONS
!$OMP SECTION
! muffin-tin exchange-correlation potential, field and energy density
call holdthd(natmtot,nthd2)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd2)
do ias=1,natmtot
  call potxcmt(tsh,ias,xctype_,rhomt_,magmt_,taumt_,exmt_,ecmt_,vxcmt_,bxcmt_, &
   wxcmt_)
end do
!$OMP END PARALLEL DO
call freethd(nthd2)
!$OMP SECTION
! interstitial exchange-correlation potential, field and energy density
call potxcir(xctype_,rhoir_,magir_,tauir_,exir_,ecir_,vxcir_,bxcir_,wxcir_)
!$OMP END SECTIONS
! symmetrise the exchange-correlation potential and magnetic field
if (tsh) then
!$OMP SECTIONS
!$OMP SECTION
  call symrf(nrmt,nrmti,npmt,ngridg,ngtot,ngvec,igfft,npmtmax,vxcmt_,vxcir_)
!$OMP SECTION
  if (spinpol) call symrvf(.true.,ncmag,nrmt,nrmti,npmt,ngridg,ngtot,ngvec, &
   igfft,npmtmax,bxcmt_,ngtot,bxcir_)
!$OMP END SECTIONS
end if
!$OMP END PARALLEL
call freethd(nthd1)
end subroutine
!EOC

