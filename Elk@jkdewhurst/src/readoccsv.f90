
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readoccsv
use modmain
implicit none
! local variables
integer ik
! get the occupation numbers from file for all k-points
do ik=1,nkpt
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
end subroutine

