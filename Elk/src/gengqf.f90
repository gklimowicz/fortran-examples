
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengqf(ng,vqpc,vgqc,gqc,jlgqr,ylmgq,sfacgq)
use modmain
implicit none
! arguments
integer, intent(in) :: ng
real(8), intent(in) :: vqpc(3)
real(8), intent(out) :: vgqc(3,ng),gqc(ng),jlgqr(njcmax,nspecies,ng)
complex(8), intent(out) :: ylmgq(lmmaxo,ng),sfacgq(ng,natmtot)
! local variables
integer ig
do ig=1,ng
! G+q-vector in Cartesian coordinates
  vgqc(:,ig)=vgc(:,ig)+vqpc(:)
! G+q-vector length
  gqc(ig)=sqrt(vgqc(1,ig)**2+vgqc(2,ig)**2+vgqc(3,ig)**2)
! spherical harmonics for G+q-vectors
  call genylmv(lmaxo,vgqc(:,ig),ylmgq(:,ig))
end do
! generate the spherical Bessel functions
call genjlgpr(ng,gqc,jlgqr)
! structure factors for G+q-vectors
call gensfacgp(ng,vgqc,ng,sfacgq)
end subroutine

