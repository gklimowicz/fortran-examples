
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine proj2d(np,fp)
use modmain
implicit none
! arguments
integer, intent(in) :: np
real(8), intent(inout) :: fp(np,3)
! local variables
integer i
real(8) vl1(3),vl2(3),t1,t2,t3
real(8) vc1(3),vc2(3),vc3(3),vc4(3)
! determine the 2D plotting plane vectors in Cartesian coordinates
vl1(:)=vclp2d(:,1)-vclp2d(:,0)
vl2(:)=vclp2d(:,2)-vclp2d(:,0)
call r3mv(avec,vl1,vc1)
call r3mv(avec,vl2,vc2)
! the z axis is orthogonal to the plotting plane vectors
call r3cross(vc1,vc2,vc3)
t1=sqrt(vc1(1)**2+vc1(2)**2+vc1(3)**2)
t2=sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
t3=sqrt(vc3(1)**2+vc3(2)**2+vc3(3)**2)
if ((t1.lt.epslat).or.(t2.lt.epslat).or.(t3.lt.epslat)) then
  write(*,*)
  write(*,'("Error(proj2d): degenerate 2D plotting directions")')
  write(*,*)
  stop
end if
! normalise the x and z axes
vc1(:)=vc1(:)/t1
vc3(:)=vc3(:)/t3
! create new y axis orthogonal to x and z axes
call r3cross(vc3,vc1,vc2)
t1=sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
vc2(:)=vc2(:)/t1
! project the vector function onto the orthogonalised plotting plane axes
do i=1,np
  vc4(:)=fp(i,:)
  fp(i,1)=dot_product(vc4(:),vc1(:))
  fp(i,2)=dot_product(vc4(:),vc2(:))
  fp(i,3)=dot_product(vc4(:),vc3(:))
end do
end subroutine

