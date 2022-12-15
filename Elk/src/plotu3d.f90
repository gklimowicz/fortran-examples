
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine plotu3d(fnum,nf,zfmt,zfir)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum,nf
complex(8), intent(in) :: zfmt(npcmtmax,natmtot,nf,nfqrz)
complex(8), intent(in) :: zfir(ngtot,nf,nfqrz)
! local variables
integer np,jf,ip
real(8) v1(3)
! allocatable arrays
real(8), allocatable :: vpl(:,:),fp(:,:)
if ((nf.lt.1).or.(nf.gt.4)) then
  write(*,*)
  write(*,'("Error(plotu3d): invalid number of functions : ",I8)') nf
  write(*,*)
  stop
end if
! total number of plot points
np=np3d(1)*np3d(2)*np3d(3)
! allocate local arrays
allocate(vpl(3,np),fp(np,nf))
! generate the 3D plotting points
call plotpt3d(vpl)
! evaluate the functions at the grid points
call plotulr(np,vpl,nf,zfmt,zfir,fp)
! write functions to file
write(fnum,'(3I6," : grid size")') np3d(:)
do ip=1,np
  call r3mv(avec,vpl(:,ip),v1)
  write(fnum,'(7G18.10)') v1(:),(fp(ip,jf),jf=1,nf)
end do
deallocate(vpl,fp)
end subroutine

