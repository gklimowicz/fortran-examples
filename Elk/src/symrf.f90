
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrf
! !INTERFACE:
subroutine symrf(nr,nri,np,ngdg,ngt,ngv,igf,ld,rfmt,rfir)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   nr   : number of radial points for each species (in,integer(nspecies))
!   nri  : number of radial points on the inner part (in,integer(nspecies))
!   np   : total number of points in each muffin-tin (in,integer(nspecies))
!   ngdg : G-vector grid sizes (in,integer(3))
!   ngt  : total number of G-vectors (in,integer)
!   ngv  : number of G-vectors within cut-off (in,integer)
!   igf  : map from G-vector index to FFT array (in,integer(ngv))
!   ld   : leading dimension (in,integer)
!   rfmt : real muffin-tin function (inout,real(ld,natmtot))
!   rfir : real intersitial function (inout,real(ngtot))
! !DESCRIPTION:
!   Symmetrises a real scalar function defined over the entire unit cell using
!   the full set of crystal symmetries. In the muffin-tin of a particular atom
!   the spherical harmonic coefficients of every equivlent atom are rotated and
!   averaged. The interstitial part of the function is first Fourier transformed
!   to $G$-space, and then averaged over each symmetry by rotating the Fourier
!   coefficients and multiplying them by a phase factor corresponding to the
!   symmetry translation.
!
! !REVISION HISTORY:
!   Created May 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
integer, intent(in) :: ngdg(3),ngt,ngv,igf(ngv),ld
real(8), intent(inout) :: rfmt(ld,natmtot),rfir(ngt)
! local variables
integer nthd
call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
call symrfmt(nr,nri,np,ld,rfmt)
!$OMP SECTION
call symrfir(ngdg,ngt,ngv,igf,rfir)
!$OMP END PARALLEL SECTIONS
call freethd(nthd)
end subroutine
!EOC

