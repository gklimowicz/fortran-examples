
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: k_tfvw1
! !INTERFACE:
elemental subroutine k_tfvw1(rho,grho2,dtdr,dtdgr2)
! !INPUT/OUTPUT PARAMETERS:
!   rho    : spin-unpolarised charge density (in,real)
!   grho2  : |grad rho|^2 (in,real)
!   dtdr   : dtau/drho (out,real)
!   dtdgr2 : dtau/d(|grad rho|^2) (out,real)
! !DESCRIPTION:
!   Calculates the derivatives $\partial\tau/\partial\rho$ and
!   $\partial\tau/\partial|\nabla\rho|^2$ of the gradient expansion of
!   the kinetic energy density $\tau$. This includes the Thomas-Fermi and
!   von Weizsacker terms:
!   $$ \tau=\frac{3}{10}(3\pi^2)^{2/3}\rho^{5/3}
!    +\frac{1}{72}\frac{|\nabla\rho|^2}{\rho}. $$
!
! !REVISION HISTORY:
!   Created December 2021 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rho,grho2
real(8), intent(out) :: dtdr,dtdgr2
! local variables
real(8), parameter :: pi=3.1415926535897932385d0
! Thomas-Fermi coefficient
real(8), parameter :: ctf=(3.d0/10.d0)*(3.d0*pi**2)**(2.d0/3.d0)
! von Weizsacker coefficient
real(8), parameter :: cvw=1.d0/72.d0
real(8) ri,t1,t2
if ((rho.lt.1.d-20).or.(grho2.lt.0.d0)) then
  dtdr=0.d0
  dtdgr2=0.d0
  return
end if
ri=1.d0/rho
t1=ctf*(5.d0/3.d0)*rho**(2.d0/3.d0)
t2=cvw*ri
dtdr=t1-t2*grho2*ri
dtdgr2=t2
end subroutine
!EOC

