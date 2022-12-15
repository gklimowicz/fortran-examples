
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gensfacgp
! !INTERFACE:
pure subroutine gensfacgp(ngp,vgpc,ld,sfacgp)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,*))
!   ld     : leading dimension (in,integer)
!   sfacgp : structure factors of G+p-vectors (out,complex(ld,natmtot))
! !DESCRIPTION:
!   Generates the atomic structure factors for a set of ${\bf G+p}$-vectors:
!   $$ S_{\alpha}({\bf G+p})=\exp(i({\bf G+p})\cdot{\bf r}_{\alpha}), $$
!   where ${\bf r}_{\alpha}$ is the position of atom $\alpha$.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: vgpc(3,ngp)
integer, intent(in) :: ld
complex(8), intent(out) :: sfacgp(ld,natmtot)
! local variables
integer is,ia,ias,igp
real(8) v1,v2,v3,t1
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  v1=atposc(1,ia,is); v2=atposc(2,ia,is); v3=atposc(3,ia,is)
  do igp=1,ngp
    t1=vgpc(1,igp)*v1+vgpc(2,igp)*v2+vgpc(3,igp)*v3
    sfacgp(igp,ias)=cmplx(cos(t1),sin(t1),8)
  end do
end do
end subroutine
!EOC

