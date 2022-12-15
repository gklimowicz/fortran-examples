
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_2a
! !INTERFACE:
subroutine ggamt_2a(tsh,is,np,rho,g2rho,gvrho,grho2)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggamt\_sp\_2a}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: is,np
real(8), intent(in) :: rho(np)
real(8), intent(out) :: g2rho(np),gvrho(np,3),grho2(np)
! local variables
integer nr,nri,i
! allocatable arrays
real(8), allocatable :: grfmt(:,:),rfmt1(:),rfmt2(:)
allocate(grfmt(np,3),rfmt2(np))
nr=nrmt(is)
nri=nrmti(is)
! compute grad^2 rho in spherical coordinates
if (tsh) then
  call grad2rfmt(nr,nri,rlmt(:,-1,is),rlmt(:,-2,is),wcrmt(:,:,is),rho,rfmt2)
else
  allocate(rfmt1(np))
  call rfsht(nr,nri,rho,rfmt1)
  call grad2rfmt(nr,nri,rlmt(:,-1,is),rlmt(:,-2,is),wcrmt(:,:,is),rfmt1,rfmt2)
end if
call rbsht(nr,nri,rfmt2,g2rho)
! compute grad rho in spherical coordinates
if (tsh) then
  call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rho,np,grfmt)
else
  call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rfmt1,np,grfmt)
end if
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),gvrho(:,i))
end do
! (grad rho)^2
grho2(1:np)=gvrho(1:np,1)**2+gvrho(1:np,2)**2+gvrho(1:np,3)**2
deallocate(rfmt2,grfmt)
if (.not.tsh) deallocate(rfmt1)
end subroutine
!EOC

