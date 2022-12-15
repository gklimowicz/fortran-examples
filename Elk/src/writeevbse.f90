
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeevbse
use modmain
implicit none
! local variables
integer ik,a
integer ios,nmbse_
! allocatable arrays
complex(8), allocatable :: w(:)
! initialise global variables
call init0
call init1
! read Fermi energy from a file
call readfermi
! get the eigenvalues from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
end do
! generate the BSE state index arrays
call genidxbse
! allocate global BSE arrays
if (allocated(evalbse)) deallocate(evalbse)
allocate(evalbse(nmbse))
if (allocated(hmlbse)) deallocate(hmlbse)
allocate(hmlbse(nmbse,nmbse))
! read in BSE Hamiltonian matrix
open(140,file='HMLBSE.OUT',form='UNFORMATTED',action='READ',status='OLD', &
 iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(writeevbse): error opening HMLBSE.OUT")')
  write(*,*)
  stop
end if
read(140) nmbse_
if (nmbse.ne.nmbse_) then
  write(*,*)
  write(*,'("Error(writeevbse): differing nmbse")')
  write(*,'(" current    : ",I6)') nmbse
  write(*,'(" HMLBSE.OUT : ",I6)') nmbse_
  write(*,*)
  stop
end if
read(140) hmlbse
close(140)
write(*,*)
write(*,'("Info(writeevbse): diagonalising the BSE Hamiltonian matrix")')
if (bsefull) then
! full non-Hermitian matrix
  allocate(w(nmbse))
  call eveqnzg(nmbse,nmbse,hmlbse,w)
  evalbse(:)=dble(w(:))
else
! Hermitian block only
  call eveqnzh(nmbse,nmbse,hmlbse,evalbse)
end if
! write the BSE eigenvectors and eigenvalues to file
open(140,file='EVBSE.OUT',form='UNFORMATTED',action='WRITE')
write(140) nmbse
write(140) evalbse
write(140) hmlbse
close(140)
! write the BSE eigenvalues to file
open(50,file='EIGVAL_BSE.OUT',form='FORMATTED',action='WRITE')
write(50,'(I6," : nmbse")') nmbse
if (bsefull) then
  do a=1,nmbse
    write(50,'(I6,2G18.10)') a,dble(w(a)),aimag(w(a))
  end do
  deallocate(w)
else
  do a=1,nmbse
    write(50,'(I6,G18.10)') a,evalbse(a)
  end do
end if
close(50)
write(*,*)
write(*,'("Info(writeevbse):")')
write(*,'(" BSE eigenvectors and eigenvalues written to EVBSE.OUT")')
write(*,'(" BSE eigenvalues written to EIGVAL_BSE.OUT")')
! deallocate global BSE arrays
deallocate(evalbse,hmlbse)
end subroutine

