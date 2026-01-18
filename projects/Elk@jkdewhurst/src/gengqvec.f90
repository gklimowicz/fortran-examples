
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengqvec(iq)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq
! local variables
integer is,ig
do ig=1,ngtot
! G+q-vector in Cartesian coordinates
  vgqc(1:3,ig)=vgc(1:3,ig)+vqc(1:3,iq)
! G+q-vector length
  gqc(ig)=sqrt(vgqc(1,ig)**2+vgqc(2,ig)**2+vgqc(3,ig)**2)
end do
! spherical harmonics for G+q-vectors
do ig=1,ngvec
  call genylmv(lmaxo,vgqc(:,ig),ylmgq(:,ig))
end do
! compute the spherical Bessel functions j_l(|G+q|R_mt)
call genjlgprmt(lnpsd,ngvec,gqc,ngvec,jlgqrmt)
! structure factors for G+q
call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! generate the smooth step function form factors for G+q
do is=1,nspecies
  call genffacgp(is,gqc,ffacgq(:,is))
end do
end subroutine

