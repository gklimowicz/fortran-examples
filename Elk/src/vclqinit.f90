
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vclqinit
use modmain
use modulr
implicit none
! zero the external Coulomb potential in Q-space
vclq(:)=0.d0
if (trdvclr) then
! read the external Coulomb potential from file if required
  call readvclr
else
! determine the external Coulomb potential from the constant electric field
  call potefieldu
! write the external Coulomb potential to file
  call writevclr
end if
end subroutine

