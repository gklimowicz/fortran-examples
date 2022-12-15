
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genvmatmt
 !INTERFACE:
subroutine genvmatmt
! !USES:
use modmain
use moddftu
! !DESCRIPTION:
!   Calculate the DFT+$U$ potential matrix to be used in the second-variational
!   step and/or the potential matrix used for fixed tensor moment calculations.
!   See {\it Phys. Rev. B} {\bf 52}, 5467 (1995) and {\it Phys. Rev. B}
!   {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created November 2007 (FC,FB,LN,JKD)
!   Fixed bug for dftu=3, January 2021 (JKD)
!   Cleaned up and removed options, September 2021 (JKD)
!EOP
!BOC
implicit none
! zero the non-local muffin-tin potential for each atom
vmatmt(:,:,:,:,:)=0.d0
! add the DFT+U potential and calculate the corresponding energies
if (dftu.ne.0) call vmatmtdu
! add the fixed tensor moment potential if required
if (ftmtype.gt.0) call vmatmtftm
! symmetrise the potential
call symdmat(lmaxdm,lmmaxdm,vmatmt)
end subroutine
!EOC

