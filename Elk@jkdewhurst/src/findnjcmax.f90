
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findnjcmax
use modmain
implicit none
! local variables
integer is,n
! find the maximum size of the spherical Bessel function array over all species
njcmax=1
do is=1,nspecies
  n=(lmaxi+1)*nrcmti(is)+(lmaxo+1)*(nrcmt(is)-nrcmti(is))
  if (n.gt.njcmax) njcmax=n
end do
end subroutine

