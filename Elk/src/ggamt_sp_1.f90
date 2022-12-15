
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_sp_1
! !INTERFACE:
subroutine ggamt_sp_1(is,np,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   is    : species number (in,integer)
!   np    : number of muffin-tin points (in,integer)
!   rhoup : spin-up density in spherical coordinates (in,real(np))
!   rhodn : spin-down density (in,real(np))
!   grho  : |grad rho| (out,real(np))
!   gup   : |grad rhoup| (out,real(np))
!   gdn   : |grad rhodn| (out,real(np))
!   g2up  : grad^2 rhoup (out,real(np))
!   g2dn  : grad^2 rhodn (out,real(np))
!   g3rho : (grad rho).(grad |grad rho|) (out,real(np))
!   g3up  : (grad rhoup).(grad |grad rhoup|) (out,real(np))
!   g3dn  : (grad rhodn).(grad |grad rhodn|) (out,real(np))
! !DESCRIPTION:
!   Computes $|\nabla\rho|$, $|\nabla\rho^{\uparrow}|$,
!   $|\nabla\rho^{\downarrow}|$, $\nabla^2\rho^{\uparrow}$,
!   $\nabla^2\rho^{\downarrow}$, $\nabla\rho\cdot(\nabla|\nabla\rho|)$,
!   $\nabla\rho^{\uparrow}\cdot(\nabla|\nabla\rho^{\uparrow}|)$ and
!   $\nabla\rho^{\downarrow}\cdot(\nabla|\nabla\rho^{\downarrow}|)$
!   for a muffin-tin charge density, as required by the generalised gradient
!   approximation functionals of type 1 for spin-polarised densities. The input
!   densities and output gradients are in terms of spherical coordinates. See
!   routines {\tt potxc} and {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created April 2004 (JKD)
!   Simplified and improved, October 2009 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is,np
real(8), intent(in) :: rhoup(np),rhodn(np)
real(8), intent(out) :: grho(np),gup(np),gdn(np)
real(8), intent(out) :: g2up(np),g2dn(np)
real(8), intent(out) :: g3rho(np),g3up(np),g3dn(np)
! local variables
integer nr,nri,i
! allocatable arrays
real(8), allocatable :: gvup(:,:),gvdn(:,:),grfmt(:,:)
real(8), allocatable :: rfmt1(:),rfmt2(:)
allocate(grfmt(np,3),gvup(np,3),gvdn(np,3))
allocate(rfmt1(np),rfmt2(np))
nr=nrmt(is)
nri=nrmti(is)
!----------------!
!     rho up     !
!----------------!
! convert rhoup to spherical harmonics
call rfsht(nr,nri,rhoup,rfmt1)
! grad rhoup in spherical coordinates
call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rfmt1,np,grfmt)
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),gvup(:,i))
end do
! |grad rhoup|
gup(1:np)=sqrt(gvup(1:np,1)**2+gvup(1:np,2)**2+gvup(1:np,3)**2)
! grad^2 rhoup in spherical coordinates
call grad2rfmt(nr,nri,rlmt(:,-1,is),rlmt(:,-2,is),wcrmt(:,:,is),rfmt1,rfmt2)
call rbsht(nr,nri,rfmt2,g2up)
! (grad rhoup).(grad |grad rhoup|)
call rfsht(nr,nri,gup,rfmt1)
call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rfmt1,np,grfmt)
g3up(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt1)
  g3up(1:np)=g3up(1:np)+gvup(1:np,i)*rfmt1(1:np)
end do
!------------------!
!     rho down     !
!------------------!
! convert rhodn to spherical harmonics
call rfsht(nr,nri,rhodn,rfmt1)
! grad rhodn in spherical coordinates
call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rfmt1,np,grfmt)
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),gvdn(:,i))
end do
gdn(1:np)=sqrt(gvdn(1:np,1)**2+gvdn(1:np,2)**2+gvdn(1:np,3)**2)
! grad^2 rhodn in spherical coordinates
call grad2rfmt(nr,nri,rlmt(:,-1,is),rlmt(:,-2,is),wcrmt(:,:,is),rfmt1,rfmt2)
call rbsht(nr,nri,rfmt2,g2dn)
! (grad rhodn).(grad |grad rhodn|)
call rfsht(nr,nri,gdn,rfmt1)
call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rfmt1,np,grfmt)
g3dn(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt1)
  g3dn(1:np)=g3dn(1:np)+gvdn(1:np,i)*rfmt1(1:np)
end do
!-------------!
!     rho     !
!-------------!
! |grad rho|
grho(1:np)=sqrt((gvup(1:np,1)+gvdn(1:np,1))**2 &
               +(gvup(1:np,2)+gvdn(1:np,2))**2 &
               +(gvup(1:np,3)+gvdn(1:np,3))**2)
! (grad rho).(grad |grad rho|)
call rfsht(nr,nri,grho,rfmt1)
call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rfmt1,np,grfmt)
g3rho(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt1)
  g3rho(1:np)=g3rho(1:np)+(gvup(1:np,i)+gvdn(1:np,i))*rfmt1(1:np)
end do
deallocate(rfmt1,rfmt2,grfmt,gvup,gvdn)
end subroutine
!EOC

