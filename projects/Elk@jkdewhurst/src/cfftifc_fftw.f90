
! Copyright (C) 2022 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine cfftifc(nd,n,sgn,c)
implicit none
! arguments
integer, intent(in) :: nd,n(nd),sgn
complex(4), intent(inout) :: c(*)
! local variables
integer, parameter :: FFTW_ESTIMATE=64
integer p
integer(8) plan
real(4) t1
! interface to FFTW version 3
!$OMP CRITICAL(cfftifc_)
call sfftw_plan_dft(plan,nd,n,c,c,sgn,FFTW_ESTIMATE)
!$OMP END CRITICAL(cfftifc_)
call sfftw_execute(plan)
!$OMP CRITICAL(cfftifc_)
call sfftw_destroy_plan(plan)
!$OMP END CRITICAL(cfftifc_)
if (sgn.eq.-1) then
  p=product(n(:))
  t1=1.e0/real(p)
  call csscal(p,t1,c,1)
end if
end subroutine

subroutine rcfftifc(nd,n,sgn,r,c)
implicit none
! arguments
integer, intent(in) :: nd,n(nd),sgn
real(4), intent(inout) :: r(*)
complex(4), intent(inout) :: c(*)
! local variables
integer, parameter :: FFTW_ESTIMATE=64
integer p
integer(8) plan
real(4) t1
!$OMP CRITICAL(rcfftifc_)
if (sgn.eq.-1) then
  call sfftw_plan_dft_r2c(plan,nd,n,r,c,FFTW_ESTIMATE)
else
  call sfftw_plan_dft_c2r(plan,nd,n,c,r,FFTW_ESTIMATE)
end if
!$OMP END CRITICAL(rcfftifc_)
call sfftw_execute(plan)
!$OMP CRITICAL(rcfftifc_)
call sfftw_destroy_plan(plan)
!$OMP END CRITICAL(rcfftifc_)
if (sgn.eq.-1) then
  p=product(n(:))
  t1=1.e0/real(p)
  p=p/n(1)
  p=p*(n(1)/2+1)
  call csscal(p,t1,c,1)
end if
end subroutine

