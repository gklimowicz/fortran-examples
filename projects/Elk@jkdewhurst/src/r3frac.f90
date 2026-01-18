
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3frac
! !INTERFACE:
pure subroutine r3frac(eps,v)
! !INPUT/OUTPUT PARAMETERS:
!   eps : zero component tolerance (in,real)
!   v   : input vector (inout,real(3))
! !DESCRIPTION:
!   Finds the fractional part of each component of a real 3-vector using the
!   function ${\rm frac}\,(x)=x-\lfloor x\rfloor$. A component is taken to be
!   zero if it lies within the intervals $[0,\epsilon)$ or $(1-\epsilon,1]$.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!   Removed iv, September 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(inout) :: v(3)
! local variables
real(8) t1
t1=1.d0-eps
v(1)=v(1)-int(v(1))
if (v(1).lt.0.d0) v(1)=v(1)+1.d0
if ((v(1).lt.eps).or.(v(1).gt.t1)) v(1)=0.d0
v(2)=v(2)-int(v(2))
if (v(2).lt.0.d0) v(2)=v(2)+1.d0
if ((v(2).lt.eps).or.(v(2).gt.t1)) v(2)=0.d0
v(3)=v(3)-int(v(3))
if (v(3).lt.0.d0) v(3)=v(3)+1.d0
if ((v(3).lt.eps).or.(v(3).gt.t1)) v(3)=0.d0
end subroutine
!EOC

