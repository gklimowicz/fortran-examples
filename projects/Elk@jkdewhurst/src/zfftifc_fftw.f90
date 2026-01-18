
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfftifc(nd,n,sgn,z)
implicit none
! arguments
integer, intent(in) :: nd,n(nd),sgn
complex(8), intent(inout) :: z(*)
! local variables
integer, parameter :: FFTW_ESTIMATE=64
integer p
integer(8) plan
real(8) t1
! interface to FFTW version 3
!$OMP CRITICAL(zfftifc_)
call dfftw_plan_dft(plan,nd,n,z,z,sgn,FFTW_ESTIMATE)
!$OMP END CRITICAL(zfftifc_)
call dfftw_execute(plan)
!$OMP CRITICAL(zfftifc_)
call dfftw_destroy_plan(plan)
!$OMP END CRITICAL(zfftifc_)
if (sgn.eq.-1) then
  p=product(n(:))
  t1=1.d0/dble(p)
  call zdscal(p,t1,z,1)
end if
end subroutine

subroutine rzfftifc(nd,n,sgn,r,z)
implicit none
! arguments
integer, intent(in) :: nd,n(nd),sgn
real(8), intent(inout) :: r(*)
complex(8), intent(inout) :: z(*)
! local variables
integer, parameter :: FFTW_ESTIMATE=64
integer p
integer(8) plan
real(8) t1
!$OMP CRITICAL(rzfftifc_)
if (sgn.eq.-1) then
  call dfftw_plan_dft_r2c(plan,nd,n,r,z,FFTW_ESTIMATE)
else
  call dfftw_plan_dft_c2r(plan,nd,n,z,r,FFTW_ESTIMATE)
end if
!$OMP END CRITICAL(rzfftifc_)
call dfftw_execute(plan)
!$OMP CRITICAL(rzfftifc_)
call dfftw_destroy_plan(plan)
!$OMP END CRITICAL(rzfftifc_)
if (sgn.eq.-1) then
  p=product(n(:))
  t1=1.d0/dble(p)
  p=p/n(1)
  p=p*(n(1)/2+1)
  call zdscal(p,t1,z,1)
end if
end subroutine

