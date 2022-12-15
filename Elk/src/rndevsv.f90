
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rndevsv(rndm,evecsv)
use modmain
use modrandom
implicit none
! arguments
real(8), intent(in) :: rndm
complex(8), intent(inout) :: evecsv(nstsv,nstsv)
! local variables
integer ist,jst
real(8) a,b
if (abs(rndm).lt.1.d-8) return
! add complex random numbers to each eigenvector
do ist=1,nstsv
  do jst=1,nstsv
    a=rndm*(randomu()-0.5d0)
    b=rndm*(randomu()-0.5d0)
    evecsv(ist,jst)=evecsv(ist,jst)+cmplx(a,b,8)
  end do
end do
! orthonormalise the eigenvectors
call unitary(nstsv,evecsv)
end subroutine

