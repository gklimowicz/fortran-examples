
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writegclq
! !INTERFACE:
subroutine writegclq
! !USES:
use modmain
! !DESCRIPTION:
!   Outputs the volume-averaged integral of $4\pi/q^2$ in the small
!   parallelepiped around each discrete $q$-point to the file {\tt GCLQ.OUT}.
!   These represent the regularised Coulomb Green's function in reciprocal
!   space for small $q$. See the routine gengclq.
!
! !REVISION HISTORY:
!   Created June 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer iq
real(8) t1
open(50,file='GCLQ'//trim(filext),form='FORMATTED',action='WRITE')
write(50,'(I6," : nqpt; q-point, vql, gclq, 4π/q² below")') nqpt
do iq=1,nqpt
  t1=vqc(1,iq)**2+vqc(2,iq)**2+vqc(3,iq)**2
  if (t1.gt.1.d-12) then
    t1=fourpi/t1
  else
    t1=0.d0
  end if
  write(50,'(I6,5G18.10)') iq,vql(:,iq),gclq(iq),t1
end do
close(50)
end subroutine
!EOC

