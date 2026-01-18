
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findkpt(vpl,isym,ik)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
integer, intent(out) :: isym,ik
! local variables
integer i1,i2,i3,lspl
real(8) v1(3),v2(3),t1
v1(:)=vpl(:)-vkloff(:)/dble(ngridk(:))
i1=modulo(nint(v1(1)*ngridk(1)),ngridk(1))
i2=modulo(nint(v1(2)*ngridk(2)),ngridk(2))
i3=modulo(nint(v1(3)*ngridk(3)),ngridk(3))
ik=ivkik(i1,i2,i3)
v1(:)=vkl(:,ik)
! find the symmetry which rotates vkl to vpl
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
! multiply vpl by the transpose of the symmetry matrix (i.e. the inverse)
  v2(:)=symlat(1,:,lspl)*vpl(1) &
       +symlat(2,:,lspl)*vpl(2) &
       +symlat(3,:,lspl)*vpl(3)
  call r3frac(epslat,v2)
  t1=abs(v1(1)-v2(1))+abs(v1(2)-v2(2))+abs(v1(3)-v2(3))
  if (t1.lt.epslat) return
end do
write(*,*)
write(*,'("Error(findkpt): equivalent k-point not in set")')
write(*,'(" Requested k-point : ",3G18.10)') vpl
write(*,*)
stop
end subroutine

