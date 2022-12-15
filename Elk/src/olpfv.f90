
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpfv(nmatp,ngp,igpig,apwalm,o)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: nmatp,ngp,igpig(ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: o(nmatp,nmatp)
! local variables
integer is,ias,j,nthd
! zero the upper triangular part of the matrix
do j=1,nmatp
  o(1:j,j)=0.d0
end do
do ias=1,natmtot
  is=idxis(ias)
  call olpaa(tefvr,is,ngp,apwalm(:,:,:,ias),nmatp,o)
end do
call olpistl(ngp,igpig,nmatp,o)
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call olpalo(is,ias,ngp,apwalm(:,:,:,ias),nmatp,o)
  call olplolo(is,ias,ngp,nmatp,o)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

