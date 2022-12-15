
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dmtotm3
! !INTERFACE:
subroutine dmtotm3(l,k,p,r,ld,dm,wkpr)
! !INPUT/OUTPUT PARAMETERS:
!   l    : angular momentum (in,integer)
!   k    : k-index of tensor moment (in,integer)
!   p    : p-index of tensor moment (in,integer)
!   r    : r-index of tensor moment (in,integer)
!   ld   : leading dimension (in,integer)
!   dm   : density matrix (in,complex(ld,2,ld,2))
!   wkpr : 3-index spherical tensor moments (out,real(-ld:ld))
! !DESCRIPTION:
!   Determines the 3-index spherical tensor moments of a density matrix $D$ with
!   $$ w_t^{kpr}=\tr\big(\Gamma_t^{kpr}D\big). $$
!   This exploits the orthonormality of the $\Gamma_t^{kpr}$ matrices. See the
!   routines {\tt tm2todm} and {\tt tm3todm} for more details.
!
! !REVISION HISTORY:
!   Created April 2008 (F. Cricchio and L. Nordstrom)
!   Modified, January 2014 (JKD)
!   Changed to real tensor moments, December 2021 (JKD)
!EOP
!BOC
implicit none
integer, intent(in) :: l,k,p,r,ld
complex(8), intent(in) :: dm(ld,2,ld,2)
real(8), intent(out) :: wkpr(-ld:ld)
! local variables
integer n,t
! automatic arrays
real(8) w(-ld:ld)
complex(8) gamma(ld,2,ld,2)
! external functions
real(8), external :: trzhmm
if (l.lt.0) then
  write(*,*)
  write(*,'("Error(dmtotm3): l < 0 : ",I8)') l
  write(*,*)
  stop
end if
if (k.lt.0) then
  write(*,*)
  write(*,'("Error(dmtotm3): k < 0 : ",I8)') k
  write(*,*)
  stop
end if
if (k.gt.2*l) then
  write(*,*)
  write(*,'("Error(dmtotm3): k > 2*l : ",2I8)') k,2*l
  write(*,*)
  stop
end if
if ((p.lt.0).or.(p.gt.1)) then
  write(*,*)
  write(*,'("Error(dmtotm3): p should be 0 or 1 : ",I8)') p
  write(*,*)
  stop
end if
if (r.lt.abs(k-p)) then
  write(*,*)
  write(*,'("Error(dmtotm3): r < |k-p| : ",2I8)') r,abs(k-p)
  write(*,*)
  stop
end if
if (r.gt.(k+p)) then
  write(*,*)
  write(*,'("Error(dmtotm3): r > k+p : ",2I8)') r,k+p
  write(*,*)
  stop
end if
n=ld*2
wkpr(:)=0.d0
do t=-r,r
  w(:)=0.d0
  w(t)=1.d0
  call tm3todm(l,k,p,r,ld,w,gamma)
  wkpr(t)=trzhmm(n,gamma,dm)
end do
end subroutine
!EOC

