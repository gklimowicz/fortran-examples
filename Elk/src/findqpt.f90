
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findqpt(vpl,isym,iq)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
integer, intent(out) :: isym,iq
! local variables
integer i1,i2,i3,lspl
real(8) v1(3),v2(3),t1
i1=modulo(nint(vpl(1)*ngridq(1)),ngridq(1))
i2=modulo(nint(vpl(2)*ngridq(2)),ngridq(2))
i3=modulo(nint(vpl(3)*ngridq(3)),ngridq(3))
iq=ivqiq(i1,i2,i3)
v1(:)=vql(:,iq)
call r3frac(epslat,v1)
! find the symmetry which rotates vql to vpl
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
write(*,'("Error(findqpt): equivalent q-point not in set")')
write(*,'(" Requested q-point : ",3G18.10)') vpl
write(*,*)
stop
end subroutine

