
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: k_tfvw
! !INTERFACE:
subroutine k_tfvw(n,rho,grho2,dtdr,dtdgr2)
! !INPUT/OUTPUT PARAMETERS:
!   n      : number of density points (in,integer)
!   rho    : spin-unpolarised charge density (in,real(n))
!   grho2  : |grad rho|^2 (in,real(n))
!   dtdr   : dtau/drho (out,real(n))
!   dtdgr2 : dtau/d|grad rho|^2 (out,real(n))
! !DESCRIPTION:
!   Calculates the derivatives $\partial\tau/\partial\rho$ and
!   $\partial\tau/\partial|\nabla\rho|^2$ of the gradient expansion of the
!   kinetic energy density $\tau$ for a set of points. See {\tt k\_tfvw1}.
!
! !REVISION HISTORY:
!   Created December 2021 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rho(n),grho2(n)
real(8), intent(out) :: dtdr(n),dtdgr2(n)
! local variables
integer i
do i=1,n
  call k_tfvw1(rho(i),grho2(i),dtdr(i),dtdgr2(i))
end do
end subroutine
!EOC

