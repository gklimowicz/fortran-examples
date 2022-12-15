
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_1
! !INTERFACE:
subroutine ggamt_1(tsh,is,np,rho,grho,g2rho,g3rho)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggamt\_sp\_1}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: is,np
real(8), intent(in) :: rho(np)
real(8), intent(out) :: grho(np),g2rho(np),g3rho(np)
! local variables
integer nr,nri,i
! allocatable arrays
real(8), allocatable :: grfmt(:,:),gvrho(:,:),rfmt1(:),rfmt2(:)
allocate(grfmt(np,3),gvrho(np,3),rfmt2(np))
nr=nrmt(is)
nri=nrmti(is)
! |grad rho|
if (tsh) then
  call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rho,np,grfmt)
else
  allocate(rfmt1(np))
  call rfsht(nr,nri,rho,rfmt1)
  call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rfmt1,np,grfmt)
end if
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),gvrho(:,i))
end do
grho(1:np)=sqrt(gvrho(1:np,1)**2+gvrho(1:np,2)**2+gvrho(1:np,3)**2)
! grad^2 rho in spherical coordinates
if (tsh) then
  call grad2rfmt(nr,nri,rlmt(:,-1,is),rlmt(:,-2,is),wcrmt(:,:,is),rho,rfmt2)
else
  call grad2rfmt(nr,nri,rlmt(:,-1,is),rlmt(:,-2,is),wcrmt(:,:,is),rfmt1,rfmt2)
end if
call rbsht(nr,nri,rfmt2,g2rho)
! (grad rho).(grad |grad rho|)
call rfsht(nr,nri,grho,rfmt2)
call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rfmt2,np,grfmt)
g3rho(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt2)
  g3rho(1:np)=g3rho(1:np)+gvrho(1:np,i)*rfmt2(1:np)
end do
deallocate(grfmt,gvrho,rfmt2)
if (.not.tsh) deallocate(rfmt1)
end subroutine
!EOC

