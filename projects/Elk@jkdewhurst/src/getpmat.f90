
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getpmat(vpl,pmat)
use modmain
use modramdisk
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: pmat(nstsv,nstsv,3)
! local variables
logical tgs
integer isym,ik,ist,jst
integer recl,nstsv_
real(8) vkl_(3),sc(3,3),t1
complex(8) v1(3),v2(3)
! find the equivalent k-point number and symmetry which rotates vkl to vpl
call findkpt(vpl,isym,ik)
! find the record length
inquire(iolength=recl) vkl_,nstsv_,pmat
!$OMP CRITICAL(u230)
! read from RAM disk if required
if (ramdisk) then
  call getrd('PMAT.OUT',ik,tgs,v1=vkl_,n1=nstsv_,nzv=nstsv*nstsv*3,zva=pmat)
  if (tgs) goto 10
end if
open(230,file='PMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(230,rec=ik) vkl_,nstsv_,pmat
close(230)
10 continue
!$OMP END CRITICAL(u230)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getpmat): differing vectors for k-point ",I8)') ik
  write(*,'(" current  : ",3G18.10)') vkl(:,ik)
  write(*,'(" PMAT.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getpmat): differing nstsv for k-point ",I8)') ik
  write(*,'(" current  : ",I8)') nstsv
  write(*,'(" PMAT.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
! if p = k then return
t1=abs(vpl(1)-vkl(1,ik))+abs(vpl(2)-vkl(2,ik))+abs(vpl(3)-vkl(3,ik))
if (t1.lt.epslat) return
! rotate the matrix elements from the reduced to non-reduced k-point
sc(:,:)=symlatc(:,:,lsplsymc(isym))
do jst=1,nstsv
  do ist=1,nstsv
    v1(:)=pmat(ist,jst,:)
    call rz3mv(sc,v1,v2)
    pmat(ist,jst,:)=v2(:)
  end do
end do
return

contains

pure subroutine rz3mv(a,x,y)
implicit none
real(8), intent(in) :: a(3,3)
complex(8), intent(in) :: x(3)
complex(8), intent(out) :: y(3)
y(1)=a(1,1)*x(1)+a(1,2)*x(2)+a(1,3)*x(3)
y(2)=a(2,1)*x(1)+a(2,2)*x(2)+a(2,3)*x(3)
y(3)=a(3,1)*x(1)+a(3,2)*x(2)+a(3,3)*x(3)
end subroutine

end subroutine

