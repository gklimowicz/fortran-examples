
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

elemental subroutine k_tf1(rho,dtdr)
implicit none
! arguments
real(8), intent(in) :: rho
real(8), intent(out) :: dtdr
! local variables
real(8), parameter :: pi=3.1415926535897932385d0
! Thomas-Fermi coefficient
real(8), parameter :: ctf=(3.d0/10.d0)*(3.d0*pi**2)**(2.d0/3.d0)
if (rho.lt.1.d-20) then
  dtdr=0.d0
  return
end if
dtdr=ctf*(5.d0/3.d0)*rho**(2.d0/3.d0)
end subroutine

