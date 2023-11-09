
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlfv(nmatp,ngp,igpig,vgpc,apwalm,h)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: nmatp,ngp,igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: h(nmatp,nmatp)
! local variables
integer is,ias,j,nthd
! zero the upper triangular part of the matrix
do j=1,nmatp
  h(1:j,j)=0.d0
end do
do ias=1,natmtot
  is=idxis(ias)
  call hmlaa(tefvr,is,ias,ngp,apwalm(:,:,:,ias),nmatp,h)
end do
call hmlistl(ngp,igpig,vgpc,nmatp,h)
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call hmlalo(is,ias,ngp,apwalm(:,:,:,ias),nmatp,h)
  call hmllolo(is,ias,ngp,nmatp,h)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

