
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine jprplot
use modmain
use modmpi
implicit none
! local variables
integer ik
! initialise universal variables
call init0
call init1
! read density and potentials from file
call readstate
! generate the core wavefunctions and densities
call gencore
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the occupancies from file
do ik=1,nkpt
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! compute the paramagnetic current density
call genjpr
! plot the current density (master process only)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(jprplot):")')
  select case(task)
  case(371)
    open(50,file='JPR1D.OUT',form='FORMATTED')
    open(51,file='JPRLINES.OUT',form='FORMATTED')
    call plot1d(50,51,3,jrmt,jrir)
    close(50)
    write(*,'(" 1D paramagnetic current density written to JPR1D.OUT")')
    write(*,'(" vertex location lines written to JPRLINES.OUT")')
  case(372)
    open(50,file='JPR2D.OUT',form='FORMATTED')
    call plot2d(.true.,50,3,jrmt,jrir)
    close(50)
    write(*,'(" 2D paramagnetic current density written to JPR2D.OUT")')
    write(*,'(" Note that the 3D vector field has been locally projected")')
    write(*,'(" onto the 2D plotting plane axes")')
  case(373)
    open(50,file='JPR3D.OUT',form='FORMATTED')
    call plot3d(50,3,jrmt,jrir)
    close(50)
    write(*,'(" 3D paramagnetic current density written to JPR3D.OUT")')
  end select
end if
end subroutine

