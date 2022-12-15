
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine dhmlistl(ngp,ngpq,igpig,igpqig,vgpc,vgpqc,ld,dh)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ngp,ngpq
integer, intent(in) :: igpig(ngkmax),igpqig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax),vgpqc(3,ngkmax)
integer, intent(in) :: ld
complex(8), intent(inout) :: dh(ld,*)
! local variables
integer j1,j2,j3
integer ig,i,j
real(8) v1,v2,v3,t1
do j=1,ngp
  ig=igpig(j)
  j1=ivg(1,ig); j2=ivg(2,ig); j3=ivg(3,ig)
  v1=0.5d0*vgpc(1,j); v2=0.5d0*vgpc(2,j); v3=0.5d0*vgpc(3,j)
  do i=1,ngpq
    ig=igpqig(i)
    ig=ivgig(ivg(1,ig)-j1,ivg(2,ig)-j2,ivg(3,ig)-j3)
    t1=vgpqc(1,i)*v1+vgpqc(2,i)*v2+vgpqc(3,i)*v3
    dh(i,j)=dh(i,j)+dvsig(ig)+t1*dcfunig(ig)
  end do
end do
end subroutine

