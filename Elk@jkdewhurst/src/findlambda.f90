
! Copyright (C) 2009 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findlambda
! !INTERFACE:
subroutine findlambda(is,l,ufix,lambda0,lambda)
use modmpi
! !INPUT/OUTPUT PARAMETERS:
!   is      : species type (in,integer)
!   l       : angular momentum (in,integer)
!   ufix    : fixed U (in,integer)
!   lambda0 : starting value for screening length  (inout,real)
!   lambda  : screening length corresponding to fixed U (out,real)
! !DESCRIPTION:
!   Find the screening length corresponding to a fixed value of $U$ by using the
!   half-interval method in the first few steps and then the more efficient
!   secant method. For $U=0$ the code automatically sets the screening length to
!   ${\rm lambdamax}=50$. This value is enough to get $F^{(k)}\sim 10^{-3}$
!   corresponding to $U\sim 0$ (that perfectly mimics a bare DFT calculation).
!
! !REVISION HISTORY:
!   Created July 2009 (Francesco Cricchio)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is,l
real(8), intent(in) :: ufix
real(8), intent(inout) :: lambda0
real(8), intent(out) :: lambda
! local variables
! max iterations in secant algorithm
integer, parameter :: maxit=100
integer it,nit
! if ufix < umin, assume lambda=lambdamax
real(8), parameter :: umin=1.d-4
! if lambda < lambdamin, perform unscreened calculation
real(8), parameter :: lambdamin=1.d-2
! max value of lambda
! lambdamax=50 is enough to get F^(k)~1.d-3 corresponding to U~0
real(8), parameter :: lambdamax=50.d0
! initial step for lambda
real(8), parameter :: dl0=0.5d0
real(8) f,fp,lambdap,dl,tol
! external functions
real(8), external :: fyukawa,fyukawa0
! small U limit
if (ufix.lt.umin) then
  lambda=lambdamax
  if (mp_mpi) write(*,'("Info(findlambda): lambda set to lambdamax : ",&
   &G18.10)') lambdamax
  return
end if
! first perform a search of lambda with half-interval method and low accuracy
! initialize values and search upward from lambda0
lambda=lambda0
dl=dl0
fp=1.d0
tol=1.d-1
nit=0
do it=1,maxit
  if (lambda.lt.lambdamin) then
! unscreened Slater parameters
    f=fyukawa0(is,l,0)-ufix
  else
! screened Slater parameters
    f=fyukawa(is,l,0,lambda)-ufix
  end if
  if ((f*fp).lt.0) dl=-0.5d0*dl
  lambdap=lambda
  lambda=lambda+dl
  fp=f
  nit=nit+1
  if (abs(f).lt.tol) goto 10
end do
10 continue
! use the found value of lambda to continue the search with secant algorithm and
! higher accuracy
tol=1.d-8
! calculate F^(0)-ufix at lambdap value
if (lambdap.lt.lambdamin) then
! unscreened Slater parameters
  fp=fyukawa0(is,l,0)-ufix
else
! screened Slater parameters
  fp=fyukawa(is,l,0,lambdap)-ufix
end if
! start secant algorithm
do it=1,maxit
! calculate F^(0)-ufix
  if (lambda.lt.lambdamin) then
! unscreened Slater parameters
    f=fyukawa0(is,l,0)-ufix
  else
! screened Slater parameters
    f=fyukawa(is,l,0,lambda)-ufix
  end if
! if requested tolerance has been reached exit the loop
  if (abs(f).lt.tol) goto 20
! update lambda with secant algorithm and roll values
  dl=-f*((lambda-lambdap)/(f-fp))
  lambdap=lambda
  lambda=lambda+dl
  fp=f
  nit=nit+1
end do
20 continue
if (nit.ge.maxit) then
  write(*,*)
  write(*,'("Error(findlambda): max number of iterations to obtain lambda &
   &reached")')
  write(*,*)
  stop
else
! update initial value for lambda for the next iteration in the SC loop
! 0.5*dl0 is enough
  lambda0=lambda-0.5d0*dl0
  if (mp_mpi) write(*,'("Info(findlambda): lambda obtained in ",I4,&
   &" iterations")') nit
end if
end subroutine
!EOC

