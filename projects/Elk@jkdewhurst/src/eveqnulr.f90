
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnulr(ik0,evecu)
use modmain
use modulr
implicit none
! arguments
integer, intent(in) :: ik0
complex(8), intent(out) :: evecu(nstulr,nstulr)
! local variables
real(8) ts0,ts1
! generate the ultra long-range Hamiltonian
call timesec(ts0)
call genhmlu(ik0,evecu)
call timesec(ts1)
!$OMP ATOMIC
timemat=timemat+ts1-ts0
! find the eigenvalues and vectors
call timesec(ts0)
call eveqnzh(nstulr,nstulr,evecu,evalu(:,ik0))
call timesec(ts1)
!$OMP ATOMIC
timesv=timesv+ts1-ts0
end subroutine

