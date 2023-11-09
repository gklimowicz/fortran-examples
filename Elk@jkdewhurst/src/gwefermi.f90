
! Copyright (C) 2018 P. Elliott, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwefermi
use modmain
use modgw
use modmpi
use modomp
implicit none
! local variables
integer, parameter :: maxit=1000
integer ik,ist,it,nthd
real(8) e0,e1,e
real(8) chg,chgk
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(gwefermi): finding the GW Fermi energy")')
end if
! find minimum and maximum eigenvalues
e0=evalsv(1,1)
e1=e0
do ik=1,nkpt
  do ist=1,nstsv
    e=evalsv(ist,ik)
    if (e.lt.e0) e0=e
    if (e.gt.e1) e1=e
  end do
end do
do it=1,maxit
  if (mp_mpi.and.(mod(it,10).eq.0)) then
    write(*,'("Info(gwefermi): done ",I4," iterations")') it
  end if
  efermi=0.5d0*(e0+e1)
  chg=0.d0
! begin parallel loop over k-points
  call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(chgk) REDUCTION(+:chg) &
!$OMP NUM_THREADS(nthd)
  do ik=1,nkpt
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
    call gwchgk(ik,chgk)
    chg=chg+chgk
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! add charge from each process and redistribute
  if (np_mpi.gt.1) then
    call mpi_allreduce(mpi_in_place,chg,1,mpi_double_precision,mpi_sum,mpicom, &
     ierror)
  end if
  if (chg.lt.chgval) then
    e0=efermi
  else
    e1=efermi
  end if
  if ((e1-e0).lt.1.d-12) return
end do
write(*,*)
write(*,'("Warning(gwefermi): could not find GW Fermi energy")')
end subroutine

