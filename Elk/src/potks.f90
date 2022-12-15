
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potks
! !INTERFACE:
subroutine potks(txc)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   txc : .true. if the exchange-correlation energy density and potentials
!         should be calculated (in,logical)
! !DESCRIPTION:
!   Computes the Kohn-Sham effective potential by adding together the Coulomb
!   and exchange-correlation potentials. Also computes the effective magnetic
!   field. See routines {\tt potcoul} and {\tt potxc}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: txc
! local variables
integer is,ias,np
real(8) ts0,ts1
call timesec(ts0)
! compute the Coulomb potential
call potcoul
! meta-GGA variables if required
if ((xcgrad.eq.3).or.(xcgrad.eq.4)) then
! generate the kinetic energy density
  call gentau
! compute the Tran-Blaha '09 constant if required
  call xc_c_tb09
end if
! compute the exchange-correlation potential and fields
if (txc) call potxc(.true.,xctype,rhomt,rhoir,magmt,magir,taumt,tauir,exmt, &
 exir,ecmt,ecir,vxcmt,vxcir,bxcmt,bxcir,wxcmt,wxcir)
! optimised effective potential exchange potential
if (xctype(1).lt.0) call oepmain
! remove the source term of the exchange-correlation magnetic field if required
if (spinpol.and.nosource) call projsbf
! effective potential from sum of Coulomb and exchange-correlation potentials
do ias=1,natmtot
  is=idxis(ias)
  np=npmt(is)
  vsmt(1:np,ias)=vclmt(1:np,ias)+vxcmt(1:np,ias)
end do
vsir(1:ngtot)=vclir(1:ngtot)+vxcir(1:ngtot)
! smooth the interstitial potential if required
call rfirsm(msmooth,vsir)
! generate the effective magnetic fields
call genbs
call timesec(ts1)
timepot=timepot+ts1-ts0
end subroutine
!EOC

