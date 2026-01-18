
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zmdet
! !INTERFACE:
complex(8) function zmdet(n,a)
! !INPUT/OUTPUT PARAMETERS:
!   n : order of matrix (in,integer)
!   a : complex square matrix (inout,complex(n,n))
! !DESCRIPTION:
!   Calculates the determinant of a complex matrix $A$ by using its $LU$
!   decomposition with partial pivoting. Let $A=PLU$ where $P$ is the
!   permutation matrix corresponding to row interchanges, then
!   \begin{align*}
!    |A|&=|P||L||U| \\
!     &=(-1)^p\prod_{i=1}^n U_{ii},
!   \end{align*}
!   where $p$ is the number of interchanges. Note that the input matrix is
!   destroyed on exit.
!
! !REVISION HISTORY:
!   Created January 2020 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(inout) :: a(n,n)
! local variables
integer i,m,info
! automatic arrays
integer ipiv(n)
! perform an LU factorisation of the input matrix
call zgetrf(n,n,a,n,ipiv,info)
! multiply diagonal elements of U together
zmdet=a(1,1)
do i=2,n
  zmdet=zmdet*a(i,i)
end do
! determine the sign from the number of row interchanges
m=1
do i=1,n
  if (ipiv(i).ne.i) m=-m
end do
if (m.eq.-1) zmdet=-zmdet
end function
!EOC

