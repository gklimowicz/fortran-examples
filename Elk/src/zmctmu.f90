
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zmctmu(tcr,l,n,a,b,ld,c)
use modomp
implicit none
! arguments
logical, intent(in) :: tcr
integer, intent(in) :: l,n
complex(8), intent(in) :: a(l,n),b(l,n)
integer, intent(in) :: ld
complex(8), intent(out) :: c(ld,*)
! local variables
integer l2,i,j,nthd
! external functions
real(8), external :: ddot
complex(8), external :: zdotc
l2=l*2
call holdthd(n,nthd)
if (tcr) then
! matrix c is real valued
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i) &
!$OMP NUM_THREADS(nthd)
  do j=1,n
    do i=1,j
      c(i,j)=c(i,j)+ddot(l2,a(:,i),1,b(:,j),1)
    end do
  end do
!$OMP END PARALLEL DO
else
! matrix c is complex valued
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i) &
!$OMP NUM_THREADS(nthd)
  do j=1,n
    do i=1,j-1
      c(i,j)=c(i,j)+zdotc(l,a(:,i),1,b(:,j),1)
    end do
    c(j,j)=c(j,j)+ddot(l2,a(:,j),1,b(:,j),1)
  end do
!$OMP END PARALLEL DO
end if
call freethd(nthd)
end subroutine

