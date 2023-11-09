
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readdatposc
use modmain
implicit none
! local variables
integer is,ia,is_,ia_,ios
open(50,file='DATPOSC.OUT',form='FORMATTED',action='READ',status='OLD', &
 iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(readdatposc): error opening DATPOSC.OUT")')
  write(*,*)
  stop
end if
do is=1,nspecies
  do ia=1,natoms(is)
    read(50,*) is_,ia_,datposc(:,:,ia,is)
    if ((is.ne.is_).or.(ia.ne.ia_)) then
      write(*,*)
      write(*,'("Error(readdatposc): species or atom number mismatch")')
      write(*,'(" internal    : ",2I4)') is,ia
      write(*,'(" DATPOSC.OUT : ",2I4)') is_,ia_
      write(*,*)
      stop
    end if
  end do
end do
close(50)
end subroutine

