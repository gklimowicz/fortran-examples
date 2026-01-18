
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zminv(n,a)
use modomp
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(inout) :: a(n,n)
! local variables
integer info,nthd,nts
! automatic arrays
integer ipiv(n)
complex(8) work(n)
! enable MKL parallelism
call holdthd(maxthdmkl,nthd)
nts=mkl_set_num_threads_local(nthd)
call zgetrf(n,n,a,n,ipiv,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(zminv): unable to invert matrix")')
  write(*,'(" ZGETRF returned INFO = ",I8)') info
  write(*,*)
  stop
end if
call zgetri(n,a,n,ipiv,work,n,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(zminv): unable to invert matrix")')
  write(*,'(" ZGETRI returned INFO = ",I8)') info
  write(*,*)
  stop
end if
nts=mkl_set_num_threads_local(0)
call freethd(nthd)
end subroutine

