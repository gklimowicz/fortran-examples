
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: k_tfvw_sp
! !INTERFACE:
subroutine k_tfvw_sp(n,rhoup,rhodn,gup2,gdn2,dtdru,dtdrd,dtdgu2,dtdgd2)
! !INPUT/OUTPUT PARAMETERS:
!   n      : number of density points (in,integer)
!   rhoup  : spin-up charge density (in,real(n))
!   rhodn  : spin-down charge density (in,real(n))
!   gup2   : |grad rhoup|^2 (in,real(n))
!   gdn2   : |grad rhodn|^2 (in,real(n))
!   dtdru  : dtauup/drhoup (out,real(n))
!   dtdrd  : dtaudn/drhodn (out,real(n))
!   dtdgu2 : dtauup/d(|grad rhoup|^2) (out,real(n))
!   dtdgu2 : dtaudn/d(|grad rhodn|^2) (out,real(n))
! !DESCRIPTION:
!   Calculates the derivatives of the spin-polarised kinetic energy density
!   $\partial\tau^{\uparrow}/\partial\rho^{\uparrow}$,
!   $\partial\tau^{\downarrow}/\partial\rho^{\downarrow}$,
!   $\partial\tau^{\uparrow}/\partial|\nabla\rho^{\uparrow}|^2$ and
!   $\partial\tau^{\downarrow}/\partial|\nabla\rho^{\downarrow}|^2$.
!   This is done by noting the relation for the kinetic energy functional
!   [G. L. Oliver and J. P. Perdew, {\it Phys. Rev. A}
!   {\bf 20}, 397 (1979)]
!   $$ T[\rho^{\uparrow},\rho^{\downarrow}]=\tfrac{1}{2}T[2\rho^{\uparrow}]
!    +\tfrac{1}{2}T[2\rho^{\downarrow}] $$
!   and taking, for example,
!   $$ \tau^{\uparrow}(\rho^{\uparrow},|\nabla\rho^{\uparrow}|^2)
!   =\tfrac{1}{2}\tau(2\rho^{\uparrow},4|\nabla\rho^{\uparrow}|^2), $$
!   where the gradient expansion of the unpolarised kinetic energy density is
!   used for $\tau$. See the routines {\tt k\_tfvw1}, {\tt ggamt\_4},
!   {\tt ggair\_4}, {\tt potxcmt}, and {\tt potxcir}.
!
! !REVISION HISTORY:
!   Created December 2021 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rhoup(n),rhodn(n)
real(8), intent(in) :: gup2(n),gdn2(n)
real(8), intent(out) :: dtdru(n),dtdrd(n)
real(8), intent(out) :: dtdgu2(n),dtdgd2(n)
! local variables
integer i
do i=1,n
  call k_tfvw1(2.d0*rhoup(i),4.d0*gup2(i),dtdru(i),dtdgu2(i))
  dtdgu2(i)=2.d0*dtdgu2(i)
end do
do i=1,n
  call k_tfvw1(2.d0*rhodn(i),4.d0*gdn2(i),dtdrd(i),dtdgd2(i))
  dtdgd2(i)=2.d0*dtdgd2(i)
end do
end subroutine
!EOC

