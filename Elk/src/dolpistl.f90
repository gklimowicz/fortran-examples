
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine dolpistl(ngp,ngpq,igpig,igpqig,ld,od)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ngp,ngpq
integer, intent(in) :: igpig(ngkmax),igpqig(ngkmax)
integer, intent(in) :: ld
complex(8), intent(inout) :: od(ld,*)
! local variables
integer i1,i2,i3,j1,j2,j3
integer ig,i,j
do j=1,ngp
  ig=igpig(j)
  j1=ivg(1,ig); j2=ivg(2,ig); j3=ivg(3,ig)
  do i=1,ngpq
    ig=igpqig(i)
    i1=ivg(1,ig)-j1; i2=ivg(2,ig)-j2; i3=ivg(3,ig)-j3
    od(i,j)=od(i,j)+dcfunig(ivgig(i1,i2,i3))
  end do
end do
end subroutine

