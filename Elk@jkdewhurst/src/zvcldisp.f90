
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zvcldisp(zvclmt)
use modmain
use modtddft
implicit none
! arguments
complex(8), intent(inout) :: zvclmt(npmtmax,natmtot)
! local variables
integer is,ia,ias,np,i
real(8) t1
do is=1,nspecies
  np=npmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do i=1,3
      t1=-datposc(i,0,ia,is)
! add the gradient of the Coulomb potential of the nucleus plus static density
      zvclmt(1:np,ias)=zvclmt(1:np,ias)+t1*gvnsmt(1:np,i,ias)
    end do
  end do
end do
end subroutine

