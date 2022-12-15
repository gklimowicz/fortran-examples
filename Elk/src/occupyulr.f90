
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine occupyulr
use modmain
use modulr
implicit none
! local variables
integer, parameter :: maxit=1000
integer ik0,ik,ist,it
real(8) e0,e1,e
real(8) chg,x,t1
! external functions
real(8), external :: stheta
! find minimum and maximum eigenvalues
e0=evalu(1,1)
e1=e0
do ik0=1,nkpt0
  do ist=1,nstulr
    e=evalu(ist,ik0)
    if (e.lt.e0) e0=e
    if (e.gt.e1) e1=e
  end do
end do
if (e0.lt.e0min) then
  write(*,*)
  write(*,'("Warning(occupyulr): minimum eigenvalue less than minimum &
   &linearisation energy : ",2G18.10)') e0,e0min
  write(*,'(" for s.c. loop ",I5)') iscl
end if
t1=1.d0/swidth
! determine the Fermi energy using the bisection method
do it=1,maxit
  efermi=0.5d0*(e0+e1)
  chg=0.d0
  do ik0=1,nkpt0
! central k-point
    ik=(ik0-1)*nkpa+1
    do ist=1,nstulr
      e=evalu(ist,ik0)
      if (e.lt.e0min) then
        occulr(ist,ik0)=0.d0
      else
        x=(efermi-e)*t1
        occulr(ist,ik0)=occmax*stheta(stype,x)
        chg=chg+wkpt(ik)*occulr(ist,ik0)
      end if
    end do
  end do
  if (chg.lt.chgval) then
    e0=efermi
  else
    e1=efermi
  end if
  if ((e1-e0).lt.1.d-12) return
end do
write(*,*)
write(*,'("Warning(occupyulr): could not find Fermi energy")')
end subroutine

