
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfinp
! !INTERFACE:
complex(8) function zfinp(zfmt1,zfir1,zfmt2,zfir2)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   zfmt1 : first complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(npcmtmax,natmtot))
!   zfir1 : first complex interstitial function in real-space
!           (in,complex(ngtc))
!   zfmt2 : second complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(npcmtmax,natmtot))
!   zfir2 : second complex interstitial function in real-space
!           (in,complex(ngtc))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions over the entire unit
!   cell. The muffin-tin functions should be stored on the coarse radial grid.
!   In the interstitial region, the integrand is multiplied with the
!   characteristic function to remove the contribution from the muffin-tin. See
!   routines {\tt zfmtinp} and {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(in) :: zfmt1(npcmtmax,natmtot),zfir1(ngtc)
complex(8), intent(in) :: zfmt2(npcmtmax,natmtot),zfir2(ngtc)
! local variables
integer is,ias,nthd
! external functions
complex(8), external :: zdotc,zfmtinp
! interstitial contribution (requires that one of the functions has been
! multiplied by the characteristic function)
zfinp=zdotc(ngtc,zfir1,1,zfir2,1)
zfinp=zfinp*omega/dble(ngtc)
! muffin-tin contribution
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) REDUCTION(+:zfinp) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  zfinp=zfinp+zfmtinp(nrcmt(is),nrcmti(is),wrcmt(:,is),zfmt1(:,ias), &
   zfmt2(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end function
!EOC

