
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readjtot(jt)
use modmain
use modtddft
implicit none
! arguments
real(8), intent(out) :: jt(3,ntimes)
! local variables
integer its,ios
real(8) times_,t1
open(50,file='JTOT_TD.OUT',form='FORMATTED',action='READ',status='OLD', &
 iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(readjtot): error opening JTOT_TD.OUT")')
  write(*,*)
  stop
end if
do its=1,ntimes-1
  read(50,*) times_,jt(:,its)
  t1=abs(times(its)-times_)
  if (t1.gt.1.d-10) then
    write(*,*)
    write(*,'("Error(readjtot): time step mismatch for step number ",I8)') its
    write(*,'(" internal    : ",G18.10)') times(its)
    write(*,'(" JTOT_TD.OUT : ",G18.10)') times_
    write(*,*)
    stop
  end if
end do
close(50)
! set current at last time step
jt(:,ntimes)=jt(:,ntimes-1)
end subroutine

