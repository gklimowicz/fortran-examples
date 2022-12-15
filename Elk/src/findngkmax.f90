
! Copyright (C) 2002-2012 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findngkmax
! !INTERFACE:
pure subroutine findngkmax(nkpt,vkc,nspnfv,vqcss,ngv,vgc,gkmax,ngkmax)
! !INPUT/OUTPUT PARAMETERS:
!   nkpt   : number of k-points (in,integer)
!   vkc    : k-point vectors in Cartesian coordinates (in,real(3,nkpt))
!   nspnfv : number of first-variational spin components: 1 normal case, 2 for
!            spin-spiral case (in,integer)
!   vqcss  : spin-spiral q-vector, not referenced if nspnfv=1 (in,integer)
!   ngv    : number of G-vectors (in,integer)
!   vgc    : G-vectors in Cartesian coordinates (in,real(3,ngv))
!   gkmax  : maximum allowed |G+k| (in,real)
!   ngkmax : maximum number of G+k-vectors over all k-points (out,integer)
! !DESCRIPTION:
!   Determines the largest number of ${\bf G+k}$-vectors with length less than
!   {\tt gkmax} over all the $k$-points. This variable is used for allocating
!   arrays.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!   Modified, August 2012 (JKD)
!   Removed modmain and added arguments, September 2012 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nkpt
real(8), intent(in) :: vkc(3,nkpt)
integer, intent(in) :: nspnfv
real(8), intent(in) :: vqcss(3)
integer, intent(in) :: ngv
real(8), intent(in) :: vgc(3,ngv),gkmax
integer, intent(out) :: ngkmax
! local variables
integer ispn,ik,n,ig
real(8) v1,v2,v3,w1,w2,w3,t0,t1
t0=gkmax**2+1.d-6
ngkmax=0
do ispn=1,nspnfv
  do ik=1,nkpt
    if (nspnfv.eq.2) then
! spin-spiral case
      if (ispn.eq.1) then
        v1=vkc(1,ik)+0.5d0*vqcss(1)
        v2=vkc(2,ik)+0.5d0*vqcss(2)
        v3=vkc(3,ik)+0.5d0*vqcss(3)
      else
        v1=vkc(1,ik)-0.5d0*vqcss(1)
        v2=vkc(2,ik)-0.5d0*vqcss(2)
        v3=vkc(3,ik)-0.5d0*vqcss(3)
      end if
    else
      v1=vkc(1,ik)
      v2=vkc(2,ik)
      v3=vkc(3,ik)
    end if
    n=0
    do ig=1,ngv
      w1=vgc(1,ig)+v1
      w2=vgc(2,ig)+v2
      w3=vgc(3,ig)+v3
      t1=w1**2+w2**2+w3**2
      if (t1.lt.t0) n=n+1
    end do
    ngkmax=max(ngkmax,n)
  end do
end do
end subroutine
!EOC

