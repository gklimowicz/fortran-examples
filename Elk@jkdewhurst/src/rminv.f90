
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rminv(n,a)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(inout) :: a(n,n)
! local variables
integer info
! automatic arrays
integer ipiv(n)
real(8) work(n)
call dgetrf(n,n,a,n,ipiv,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(rminv): unable to invert matrix")')
  write(*,'(" DGETRF returned INFO = ",I8)') info
  write(*,*)
  stop
end if
call dgetri(n,a,n,ipiv,work,n,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(rminv): unable to invert matrix")')
  write(*,'(" DGETRI returned INFO = ",I8)') info
  write(*,*)
  stop
end if
end subroutine

