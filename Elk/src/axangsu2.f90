
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: axangsu2
subroutine axangsu2(v,th,su2)
! !INPUT/OUTPUT PARAMETERS:
!   v   : rotation axis vector (in,real(3))
!   th  : rotation angle (in,real)
!   su2 : SU(2) representation of rotation (out,complex(2,2))
! !DESCRIPTION:
!   Finds the complex ${\rm SU}(2)$ representation of a rotation defined by an
!   axis vector $\hat{\bf v}$ and angle $\theta$. The spinor rotation matrix is
!   given explicitly by
!   $$ R^{1/2}(\hat{\bf v},\theta)=I\cos\frac{\theta}{2}
!    -i(\hat{\bf v}\cdot\vec{\sigma})\sin\frac{\theta}{2}. $$
!
! !REVISION HISTORY:
!   Created August 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: v(3),th
complex(8), intent(out) :: su2(2,2)
! local variables
real(8) x,y,z,cs,sn,t1
x=v(1); y=v(2); z=v(3)
t1=sqrt(x**2+y**2+z**2)
if (t1.lt.1.d-8) then
  write(*,*)
  write(*,'("Error(axangsu2): zero length axis vector")')
  write(*,*)
  stop
end if
! normalise the vector
t1=1.d0/t1
x=x*t1; y=y*t1; z=z*t1
cs=cos(0.5d0*th)
sn=sin(0.5d0*th)
su2(1,1)=cmplx(cs,-z*sn,8)
su2(2,1)=cmplx(y*sn,-x*sn,8)
su2(1,2)=cmplx(-y*sn,-x*sn,8)
su2(2,2)=cmplx(cs,z*sn,8)
end subroutine
!EOC

