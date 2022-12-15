
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

module modgw

! maximum Matsubara frequency for the GW calculation
real(8) wmaxgw
! maximum number of Matsubara frequencies
integer nwgw
! integer grid intervals for Matsubara frequencies
integer intwgw(2)
! map from frequency index to FFT array
integer, allocatable :: iwfft(:)
! maximum fermionic Matsubara frequency index to be used for the GW calculation
integer nwfm
! maximum bosonic frequency index
integer nwbs
! imaginary frequencies used for the GW calculation
real(8), allocatable :: wgw(:)
! complex fermionic frequencies
complex(8), allocatable :: wfm(:)
! twdiag is .true. if the screened interaction W is taken to be diagonal
logical twdiag
! tsediag is .true. if the GW self-energy is taken to be diagonal
logical tsediag
! type of analytic continuation to be used for determining the self-energy on
! the real axis
integer actype
! number of poles used for fitting the self-energy matrix elements
integer npole
! number of complex shifts used in averaging the Pade approximant for the
! analytic continuation of the self-energy to the real axis
integer nspade

end module

