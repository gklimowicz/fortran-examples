
! Copyright (C) 2020 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getephmkq(iqp,ikp,ephmat)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iqp,ikp
complex(4), intent(out) :: ephmat(nstsv,nstsv,nbph)
! local variables
integer isym,lspl,iq,ik,iv(3)
if (iqp.le.nqpt) then
! q-point is in the reduced set
  iq=iqp
  ik=ikp
else
! q-point is not in the reduced set
  call findqpt(vql(:,iqp),isym,iq)
  lspl=lsplsymc(isym)
  call i3mtv(symlat(:,:,lspl),ivk(:,ikp),iv)
  iv(:)=modulo(iv(:),ngridk(:))
  ik=ivkiknr(iv(1),iv(2),iv(3))
end if
ephmat(:,:,:)=ephmkq(:,:,:,ik,iq)
end subroutine

