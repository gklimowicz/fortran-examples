
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readforcet
use modmain
use modtddft
implicit none
! local variables
integer ios,its,its_
integer is,ia,ias,is_,ia_
real(8) times_,t1
if (allocated(forcet)) deallocate(forcet)
allocate(forcet(3,natmtot,ntimes))
! read in the time-dependent total atomic forces
open(50,file='FORCETOT_TD.OUT',form='FORMATTED',action='READ',status='OLD', &
 iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(readforcet): error opening FORCETOT_TD.OUT")')
  write(*,*)
  stop
end if
do its=1,ntimes-1
  read(50,*) its_,times_
  if (its.ne.its_) then
    write(*,*)
    write(*,'("Error(readforcet): time step number mismatch")')
    write(*,'(" internal        : ",I8)') its
    write(*,'(" FORCETOT_TD.OUT : ",I8)') its_
    write(*,*)
    stop
  end if
  t1=abs(times(its)-times_)
  if (t1.gt.1.d-10) then
    write(*,*)
    write(*,'("Error(readforcet): time step mismatch for step number ",I8)') its
    write(*,'(" internal        : ",G18.10)') times(its)
    write(*,'(" FORCETOT_TD.OUT : ",G18.10)') times_
    stop
  end if
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      read(50,*) is_,ia_,forcet(:,ias,its)
      if ((is.ne.is_).or.(ia.ne.ia_)) then
        write(*,*)
        write(*,'("Error(readforcet): species or atom number mismatch for time &
         &step number ",I8)') its
        write(*,'(" internal        : ",2I4)') is,ia
        write(*,'(" FORCETOT_TD.OUT : ",2I4)') is_,ia_
        write(*,*)
        stop
      end if
    end do
  end do
end do
close(50)
! set force at last time step
forcet(:,:,ntimes)=forcet(:,:,ntimes-1)
end subroutine

