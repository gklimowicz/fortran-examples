
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlistl
! !INTERFACE:
pure subroutine hmlistl(ngp,igpig,vgpc,ld,h)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   ld    : leading dimension of h (in,integer)
!   h     : Hamiltonian matrix (inout,complex(*))
! !DESCRIPTION:
!   Computes the interstitial contribution to the Hamiltonian matrix for the APW
!   basis functions. The Hamiltonian is given by
!   $$ H^{\rm I}({\bf G+k,G'+k})=\frac{1}{2}({\bf G+k})\cdot({\bf G'+k})
!    \tilde{\Theta}({\bf G-G'})+V_s({\bf G-G'}), $$
!   where $V_s$ is the interstitial Kohn-Sham potential and $\tilde{\Theta}$ is
!   the characteristic function. See routine {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp,igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
integer, intent(in) :: ld
complex(8), intent(inout) :: h(ld,*)
! local variables
integer ig,j1,j2,j3,i,j
real(8) v1,v2,v3,t1
do j=1,ngp
  ig=igpig(j)
  j1=ivg(1,ig); j2=ivg(2,ig); j3=ivg(3,ig)
  v1=0.5d0*vgpc(1,j); v2=0.5d0*vgpc(2,j); v3=0.5d0*vgpc(3,j)
  do i=1,j
    ig=igpig(i)
    ig=ivgig(ivg(1,ig)-j1,ivg(2,ig)-j2,ivg(3,ig)-j3)
    t1=vgpc(1,i)*v1+vgpc(2,i)*v2+vgpc(3,i)*v3
    h(i,j)=h(i,j)+vsig(ig)+t1*cfunig(ig)
  end do
end do
end subroutine
!EOC

