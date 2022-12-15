
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: engyfdu
! !INTERFACE:
subroutine engyfdu(idu)
! !USES:
use modmain
use moddftu
use modmpi
! !INPUT/OUTPUT PARAMETERS:
!   idu : DFT+U entry (in,integer)
! !DESCRIPTION:
!   Calculates the energies of radial functions to be used to calculate the
!   Slater integrals. By convention those energies are chosen to be the ones at
!   the center of the band.
!
! !REVISION HISTORY:
!   Created April 2008 (F. Cricchio)
!EOP
!BOC
implicit none
integer, intent(in) :: idu
! local variables
integer is,ia,ja,ias,jas
integer nr,nri,nnf,l
logical fnd
! automatic arrays
logical done(natmmax)
real(8) vr(nrmtmax)
nnf=0
is=isldu(1,idu)
l=isldu(2,idu)
nr=nrmt(is)
nri=nrmti(is)
done(:)=.false.
do ia=1,natoms(is)
  if (done(ia)) cycle
  ias=idxas(ia,is)
  call rfmtlm(1,nr,nri,vsmt(:,ias),vr)
  vr(1:nr)=vr(1:nr)*y00
! find the center of the band starting from -0.5 Ha
  efdu(l,ias)=-0.5d0
  call findband(solsc,l,nrmt(is),rsp(1,is),vr,epsband,demaxbnd,efdu(l,ias),fnd)
  if (.not.fnd) nnf=nnf+1
  done(ia)=.true.
! copy to equivalent atoms
  do ja=1,natoms(is)
    if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
      jas=idxas(ja,is)
      efdu(l,jas)=efdu(l,ias)
      done(ja)=.true.
    end if
  end do
! end loops over atoms and species
end do
if (mp_mpi.and.(nnf.gt.0)) then
  write(*,*)
  write(*,'("Warning(engyfdu): could not find ",I3," energies")') nnf
end if
end subroutine
!EOC

