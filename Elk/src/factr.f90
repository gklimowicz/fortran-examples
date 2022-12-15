
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: factr
! !INTERFACE:
real(8) function factr(n,d)
! !INPUT/OUTPUT PARAMETERS:
!   n : numerator (in,integer)
!   d : denominator (in,integer)
! !DESCRIPTION:
!   Returns the ratio $n!/d!$ for $n,d\ge 0$. Performs no under- or overflow
!   checking.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n,d
! local variables
integer i
! external functions
real(8), external :: factn
if (d.eq.1) then
  factr=factn(n)
  return
end if
if ((n.lt.0).or.(d.lt.0)) then
  factr=0.d0
  return
end if
if (n.lt.d) then
  factr=dble(n+1)
  do i=n+2,d
    factr=factr*dble(i)
  end do
  factr=1.d0/factr
else if (n.eq.d) then
  factr=1.d0
else
  factr=dble(d+1)
  do i=d+2,n
    factr=factr*dble(i)
  end do
end if
end function
!EOC

