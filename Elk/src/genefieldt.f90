
! Copyright (C) 2020 Peter Elliott, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genefieldt
use modmain
use modtddft
implicit none
! local variables
integer np,it,i,j
real(8) t1
! automatic arrays
real(8) ya(4)
! external functions
real(8), external :: polynm
! determine the electric field at the current time step
t1=-1.d0/solsc
np=min(4,itimes)
it=itimes-np+1
do i=1,3
  do j=1,np
    ya(j)=afieldt(i,it+j-1)
  end do
  efieldt(i)=t1*polynm(1,np,times(it),ya,times(itimes))
end do
end subroutine

