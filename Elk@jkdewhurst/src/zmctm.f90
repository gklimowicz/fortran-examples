
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zmctm(l,m,n,a,b,ld,c)
use modomp
implicit none
! arguments
integer, intent(in) :: l,m,n
complex(8), intent(in) :: a(l,m),b(l,n)
integer, intent(in) :: ld
complex(8), intent(out) :: c(ld,*)
! local variables
integer i,j,nthd
! external functions
complex(8), external :: zdotc
call holdthd(n,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i) &
!$OMP NUM_THREADS(nthd)
do j=1,n
  do i=1,m
    c(i,j)=c(i,j)+zdotc(l,a(:,i),1,b(:,j),1)
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

