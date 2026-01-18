
! Copyright (C) 2020 Peter Elliott, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tdrestart
use modmain
use modtddft
use moddftu
implicit none
! last time step for which the eigenvectors were written to file
call readtimes
! generate the density and magnetisation
call rhomag
! compute the time-dependent Kohn-Sham potentials and magnetic fields
call potkst
! DFT+U
if (dftu.ne.0) then
  call gendmatmt
  call genvmatmt
end if
! read in the induced A-field and its time derivative
if (tafindt) call readafindt
! read in the atomic displacements and velocities
if (tdatpos) trddatpos=.true.
! begin TDDFT loop with next time step
itimes0=itimes0+1
end subroutine

