
! Copyright (C) 2008  F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genfdufr
! !INTERFACE:
subroutine genfdufr(idu)
! !USES:
use modmain
use moddftu
! !INPUT/OUTPUT PARAMETERS:
!   idu : DFT+U entry (in,integer)
! !DESCRIPTION:
!   Generates the radial functions used to calculate the Slater integrals
!   through a Yukawa potential.
!
! !REVISION HISTORY:
!   Created April 2008 from genapwfr (Francesco Cricchio)
!EOP
!BOC
implicit none
integer, intent(in) :: idu
! local variables
integer is,ia,ias
integer nr,nri,ir
integer nn,l
real(8) t1
! automatic arrays
real(8) vr(nrmtmax),fr(nrmtmax)
real(8) p0(nrmtmax),p1(nrmtmax),q0(nrmtmax),q1(nrmtmax)
! external functions
real(8), external :: splint
is=isldu(1,idu)
l=isldu(2,idu)
nr=nrmt(is)
nri=nrmti(is)
do ia=1,natoms(is)
  ias=idxas(ia,is)
  call rfmtlm(1,nr,nri,vsmt(:,ias),vr)
  vr(1:nr)=vr(1:nr)*y00
! integrate the radial Schrodinger equation
  call rschrodint(solsc,l,efdu(l,ias),nr,rlmt(:,1,is),vr,nn,p0,p1,q0,q1)
! normalise radial functions
  fr(1:nr)=p0(1:nr)**2
  t1=splint(nr,rlmt(:,1,is),fr)
  if (t1.lt.1.d-20) then
    write(*,*)
    write(*,'("Error(genfdufr): degenerate radial functions")')
    write(*,'(" for species ",I4)') is
    write(*,'(" atom ",I4)') ia
    write(*,'(" and angular momentum ",I4)') l
    write(*,*)
    stop
  end if
  t1=1.d0/sqrt(abs(t1))
  p0(1:nr)=t1*p0(1:nr)
! divide by r and store in global array
  do ir=1,nr
    fdufr(ir,l,ias)=p0(ir)/rsp(ir,is)
  end do
end do
end subroutine
!EOC

