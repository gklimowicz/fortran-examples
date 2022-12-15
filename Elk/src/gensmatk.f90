
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gensmatk(vpl,smat)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: smat(2,2,nstsv,nstsv)
! local variables
integer ist,jst,n
! allocatable arrays
complex(8), allocatable :: evecsv(:,:)
! external functions
complex(8), external :: zdotc
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(gensmatk): spin-unpolarised calculation")')
  write(*,*)
  stop
end if
allocate(evecsv(nstsv,nstsv))
! get the second-variational eigenvectors from file
call getevecsv(filext,0,vpl,evecsv)
n=nstfv+1
do jst=1,nstsv
  do ist=1,jst
    smat(1,1,ist,jst)=zdotc(nstfv,evecsv(1,ist),1,evecsv(1,jst),1)
    smat(2,1,ist,jst)=zdotc(nstfv,evecsv(n,ist),1,evecsv(1,jst),1)
    smat(1,2,ist,jst)=zdotc(nstfv,evecsv(1,ist),1,evecsv(n,jst),1)
    smat(2,2,ist,jst)=zdotc(nstfv,evecsv(n,ist),1,evecsv(n,jst),1)
  end do
end do
deallocate(evecsv)
end subroutine

