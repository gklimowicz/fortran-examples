
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnfvz(nmatp,h,o,evalfv,evecfv)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: nmatp
complex(8), intent(in) :: h(*),o(*)
real(8), intent(out) :: evalfv(nstfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv)
! local variables
integer lwork,info
integer m,i,nthd,nts
real(8) vl,vu
real(8) ts0,ts1
! automatic arrays
integer iwork(5*nmatp),ifail(nmatp)
! allocatable arrays
real(8), allocatable :: w(:),rwork(:)
complex(8), allocatable :: work(:)
call timesec(ts0)
allocate(w(nmatp),rwork(7*nmatp))
lwork=2*nmatp
allocate(work(lwork))
! enable MKL parallelism
call holdthd(maxthdmkl,nthd)
nts=mkl_set_num_threads_local(nthd)
! diagonalise the matrix
call zhegvx(1,'V','I','U',nmatp,h,nmatp,o,nmatp,vl,vu,1,nstfv,evaltol,m,w, &
 evecfv,nmatmax,work,lwork,rwork,iwork,ifail,info)
nts=mkl_set_num_threads_local(0)
call freethd(nthd)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(eveqnfvz): diagonalisation failed")')
  write(*,'(" ZHEGVX returned INFO = ",I8)') info
  if (info.gt.nmatp) then
    i=info-nmatp
    write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
    write(*,'("  is not positive definite")')
    write(*,'(" Order of overlap matrix : ",I8)') nmatp
  end if
  write(*,*)
  stop
end if
evalfv(1:nstfv)=w(1:nstfv)
deallocate(w,rwork,work)
call timesec(ts1)
!$OMP ATOMIC
timefv=timefv+ts1-ts0
end subroutine

