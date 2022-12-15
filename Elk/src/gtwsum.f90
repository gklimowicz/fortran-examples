
! Copyright (C) 2018 A. Davydov, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure complex(8) function gtwsum(gf,gl)
use modmain
use modgw
implicit none
! arguments
complex(8), intent(in) :: gf,gl
! local variables
integer iw
real(8) b
complex(8) a1,a2,z1
z1=wfm(0)
a1=0.5d0*(gf-gl)*z1
a2=0.5d0*(gf+gl)*z1**2
b=1.d0/(kboltz*tempk)
gtwsum=b*(0.5d0*a1-0.25d0*b*a2)
do iw=0,nwfm
  gtwsum=gtwsum-a2/wfm(iw)**2
end do
end function

