
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnzh(n,ld,a,w)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: n,ld
complex(8), intent(inout) :: a(ld,n)
real(8), intent(out) :: w(n)
! local variables
integer lrwork,lwork,info
integer nthd,nts
! automatic arrays
integer iwork(5*n+3)
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)
! use the divide-and-conquer LAPACK routine zheevd
lrwork=2*n**2+5*n+1
lwork=n**2+2*n
allocate(rwork(lrwork),work(lwork))
! enable MKL parallelism
call holdthd(maxthdmkl,nthd)
nts=mkl_set_num_threads_local(nthd)
call zheevd('V','U',n,a,ld,w,work,lwork,rwork,lrwork,iwork,5*n+3,info)
nts=mkl_set_num_threads_local(0)
call freethd(nthd)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(eveqnzh): diagonalisation failed")')
  write(*,'(" ZHEEVD returned INFO = ",I8)') info
  write(*,*)
  stop
end if
deallocate(rwork,work)
end subroutine

