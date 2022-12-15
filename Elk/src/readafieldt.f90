
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readafieldt
use modmain
use modtddft
implicit none
! local variables
integer ios,ntimes_,its,its_
real(8) times_,t1
! generate the time step grid
call gentimes
! read in the time-dependent vector potential
open(50,file='AFIELDT.OUT',form='FORMATTED',action='READ',status='OLD', &
 iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(readafieldt): error opening AFIELDT.OUT")')
  write(*,*)
  stop
end if
read(50,*) ntimes_
if (ntimes_.le.0) then
  write(*,*)
  write(*,'("Error(readafieldt): ntimes <= 0 : ",I8)') ntimes_
  write(*,*)
  stop
end if
ntimes=min(ntimes,ntimes_)
if (allocated(afieldt)) deallocate(afieldt)
allocate(afieldt(3,ntimes))
do its=1,ntimes
  read(50,*) its_,times_,afieldt(:,its)
  if (its.ne.its_) then
    write(*,*)
    write(*,'("Error(readafieldt): time step number mismatch")')
    write(*,'(" internal    : ",I8)') its
    write(*,'(" AFIELDT.OUT : ",I8)') its_
    write(*,*)
    stop
  end if
  t1=abs(times(its)-times_)
  if (t1.gt.1.d-10) then
    write(*,*)
    write(*,'("Error(readafieldt): time step mismatch for step number ",&
     &I8)') its
    write(*,'(" internal    : ",G18.10)') times(its)
    write(*,'(" AFIELDT.OUT : ",G18.10)') times_
    stop
  end if
end do
close(50)
end subroutine

