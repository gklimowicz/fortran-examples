
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: trzhmm
! !INTERFACE:
pure real(8) function trzhmm(n,a,b)
! !INPUT/OUTPUT PARAMETERS:
!   n  : order of matrix (in,integer)
!   a  : Hermitian matrix A (in,complex(n,n))
!   b  : Hermitian matrix B (in,complex(n,n))
! !DESCRIPTION:
!   Calculates the trace of the product of two Hermitian matrices, $\tr(AB)$.
!   Only the upper triangular parts of $A$ and $B$ are referenced.
!
! !REVISION HISTORY:
!   Created December 2021 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: a(n,n),b(n,n)
! local variables
integer i,j
real(8) sm
sm=0.d0
! off-diagonal contribution (use upper triangular part)
do j=1,n
  do i=1,j-1
    sm=sm+dble(a(i,j)*conjg(b(i,j)))
  end do
end do
sm=sm*2.d0
! diagonal contribution
do i=1,n
  sm=sm+dble(a(i,i))*dble(b(i,i))
end do
trzhmm=sm
end function
!EOC

