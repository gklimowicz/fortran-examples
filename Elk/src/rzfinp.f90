
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function rzfinp(rfmt,rfir,zfmt,zfir)
use modmain
use modomp
implicit none
! arguments
real(8), intent(in) :: rfmt(npcmtmax,natmtot),rfir(ngtot)
complex(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtot)
! local variables
integer is,ias,ir,nthd
! external functions
complex(8), external :: rzfmtinp
! interstitial contribution
rzfinp=0.d0
do ir=1,ngtot
  rzfinp=rzfinp+(cfunir(ir)*rfir(ir))*zfir(ir)
end do
rzfinp=rzfinp*(omega/dble(ngtot))
! muffin-tin contribution
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) REDUCTION(+:rzfinp) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  rzfinp=rzfinp+rzfmtinp(nrcmt(is),nrcmti(is),wrcmt(:,is),rfmt(:,ias), &
   zfmt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end function

