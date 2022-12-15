
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine plotu1d(fnum1,fnum2,nf,zfmt,zfir)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum1,fnum2,nf
complex(8), intent(in) :: zfmt(npcmtmax,natmtot,nf,nfqrz)
complex(8), intent(in) :: zfir(ngtot,nf,nfqrz)
! local variables
integer jf,ip,iv
real(8) fmin,fmax,t1
! allocatable arrays
real(8), allocatable :: fp(:,:)
if ((nf.lt.1).or.(nf.gt.4)) then
  write(*,*)
  write(*,'("Error(plotu1d): invalid number of functions : ",I8)') nf
  write(*,*)
  stop
end if
allocate(fp(npp1d,nf))
! connect the 1D plotting vertices
call plotpt1d(avec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
! evaluate function at each point
call plotulr(npp1d,vplp1d,nf,zfmt,zfir,fp)
do ip=1,npp1d
! write the point distances and function to file
  write(fnum1,'(5G18.10)') dpp1d(ip),(fp(ip,jf),jf=1,nf)
end do
! write the vertex location lines
fmin=minval(fp(:,:))
fmax=maxval(fp(:,:))
t1=0.5d0*(fmax-fmin)
fmin=fmin-t1
fmax=fmax+t1
do iv=1,nvp1d
  write(fnum2,'(2G18.10)') dvp1d(iv),fmin
  write(fnum2,'(2G18.10)') dvp1d(iv),fmax
  write(fnum2,*)
end do
deallocate(fp)
end subroutine

