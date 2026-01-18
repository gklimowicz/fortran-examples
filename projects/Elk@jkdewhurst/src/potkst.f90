
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potkst
use modmain
use modtddft
implicit none
! adiabatic approximation to the exchange-correlation potential
call potks(.true.)
end subroutine

