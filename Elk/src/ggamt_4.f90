
! Copyright (C) 2020 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ggamt_4(is,np,gvrho,dtdr,dtdgr2,wx,wc,vx,vc)
use modmain
implicit none
! arguments
integer, intent(in) :: is,np
real(8), intent(in) :: gvrho(np,3)
real(8), intent(in) :: dtdr(np),dtdgr2(np)
real(8), intent(in) :: wx(np),wc(np)
real(8), intent(inout) :: vx(np),vc(np)
! local variables
integer nr,nri,i
! automatic arrays
real(8) grfmt(np,3),rfmt1(np),rfmt2(np),rfmt3(np)
nr=nrmt(is)
nri=nrmti(is)
!------------------!
!     exchange     !
!------------------!
vx(1:np)=vx(1:np)+wx(1:np)*dtdr(1:np)
rfmt1(1:np)=wx(1:np)*dtdgr2(1:np)
do i=1,3
  rfmt2(1:np)=rfmt1(1:np)*gvrho(1:np,i)
  call rfsht(nr,nri,rfmt2,rfmt3)
  call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rfmt3,np,grfmt)
  call rbsht(nr,nri,grfmt(:,i),rfmt2)
  vx(1:np)=vx(1:np)-2.d0*rfmt2(1:np)
end do
!---------------------!
!     correlation     !
!---------------------!
vc(1:np)=vc(1:np)+wc(1:np)*dtdr(1:np)
rfmt1(1:np)=wc(1:np)*dtdgr2(1:np)
do i=1,3
  rfmt2(1:np)=rfmt1(1:np)*gvrho(1:np,i)
  call rfsht(nr,nri,rfmt2,rfmt3)
  call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rfmt3,np,grfmt)
  call rbsht(nr,nri,grfmt(:,i),rfmt2)
  vc(1:np)=vc(1:np)-2.d0*rfmt2(1:np)
end do
end subroutine

