
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genjlgprmt
! !INTERFACE:
subroutine genjlgprmt(lmax,ngp,gpc,ld,jlgprmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lmax    : angular momentum cut-off (in,integer)
!   ngp     : number of G+p-vectors (in,integer)
!   gpc     : length of G+p-vectors (in,real(ngkmax))
!   ld      : leading dimension (in,integer)
!   jlgprmt : spherical Bessel functions (out,real(0:lmax,ld,nspecies))
! !DESCRIPTION:
!   Calculates and stores the spherical Bessel functions
!   $j_l(|{\bf G}+{\bf p}|{\bf R}_{\rm MT})$ for all input ${\bf G}+{\bf p}$
!   vectors and the muffin-tin radii ${\bf R}_{\rm MT}$ of every atomic species.
!
! !REVISION HISTORY:
!   Created April 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax,ngp
real(8), intent(in) :: gpc(ngp)
integer, intent(in) :: ld
real(8), intent(out) :: jlgprmt(0:lmax,ld,nspecies)
! local variables
integer is,ig
real(8) r,t1
do is=1,nspecies
  r=rmt(is)
  do ig=1,ngp
    t1=gpc(ig)*r
    call sbessel(lmax,t1,jlgprmt(:,ig,is))
  end do
end do
end subroutine
!EOC

