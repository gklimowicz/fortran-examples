
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine plotu2d(tproj,fnum,nf,zfmt,zfir)
use modmain
implicit none
! arguments
logical, intent(in) :: tproj
integer, intent(in) :: fnum,nf
complex(8), intent(in) :: zfmt(npcmtmax,natmtot,nf,nfqrz)
complex(8), intent(in) :: zfir(ngtot,nf,nfqrz)
! local variables
integer np,jf,ip
real(8) vpnl(3)
! allocatable arrays
real(8), allocatable :: vpl(:,:),vppc(:,:),fp(:,:)
if ((nf.lt.1).or.(nf.gt.4)) then
  write(*,*)
  write(*,'("Error(plotu2d): invalid number of functions : ",I8)') nf
  write(*,*)
  stop
end if
! allocate local arrays
np=np2d(1)*np2d(2)
allocate(vpl(3,np),vppc(2,np),fp(np,nf))
! generate the 2D plotting points
call plotpt2d(avec,ainv,vpnl,vpl,vppc)
! evaluate the functions at the grid points
call plotulr(np,vpl,nf,zfmt,zfir,fp)
! project the vector function onto the 2D plotting plane if required
if (tproj.and.(nf.eq.3)) then
  call proj2d(np,fp)
end if
! write the functions to file
write(fnum,'(2I6," : grid size")') np2d(:)
do ip=1,np
  write(fnum,'(6G18.10)') vppc(1,ip),vppc(2,ip),(fp(ip,jf),jf=1,nf)
end do
deallocate(vpl,vppc,fp)
end subroutine

