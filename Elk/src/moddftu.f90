
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module moddftu
use modmain

!-----------------------------------------------------------!
!     muffin-tin density and potential matrix variables     !
!-----------------------------------------------------------!
! maximum angular momentum for muffin-tin density matrix
integer, parameter :: lmaxdm=3
integer, parameter :: lmmaxdm=(lmaxdm+1)**2
! density matrix in each muffin-tin
complex(8), allocatable :: dmatmt(:,:,:,:,:)
! tolerance for checking invariance of density matrix under symmetry operations
real(8) epsdmat
! potential matrix in each muffin-tin
complex(8), allocatable :: vmatmt(:,:,:,:,:)
! potential matrix in spherical coordinates for lmaxi and lmaxo
complex(4), allocatable :: vmatmti(:,:,:,:,:),vmatmto(:,:,:,:,:)
! tvmatmt is .true. if the potential matrices are calculated
logical tvmatmt
! tvmmt is .true. if the potential matrix for that l and atom is non-zero
logical, allocatable :: tvmmt(:,:)

!-------------------------!
!     DFT+U variables     !
!-------------------------!
! type of DFT+U to use (0 = none)
integer dftu
! input type for DFT+U calculation (1:5)
integer inpdftu
! maximum number of DFT+U entries
integer, parameter :: maxdftu=10
! number of DFT+U entries
integer ndftu
! species and angular momentum for each entry
integer isldu(2,maxdftu)
! U and J values for each entry
real(8) ujdu(2,maxdftu)
! DFT+U energy for each atom and entry
real(8), allocatable :: engyadu(:,:)
! energy from the DFT+U correction
real(8) engydu
! Slater parameters
real(8) fdu(0:2*lmaxdm,maxdftu)
! Racah parameters
real(8) edu(0:lmaxdm,maxdftu)
! screening length of Yukawa potential to calculate Slater integrals
real(8) lamdu(maxdftu)
! initial values of screening length if U is fixed
real(8) lamdu0(maxdftu)
! energies to calculate radial functions for Slater integrals
real(8), allocatable :: efdu(:,:)
! radial functions to calculate Slater integrals
real(8), allocatable :: fdufr(:,:,:)
! fixed value of U for which screening length has to be determined
real(8) udufix(maxdftu),dudufix(maxdftu)

!---------------------------------!
!     tensor moment variables     !
!---------------------------------!
! tmwrite is .true. if tensor moments are written out at every s.c. loop
logical tmwrite
! fixed tensor moment type
!  0      : none
!  1 (-1) : fixed 3-index tensor moment (or just lowering the symmetry)
integer ftmtype
! number of fixed tensor moment entries
integer ntmfix
! tensor moment indices for each entry:
!  is, ia, l for the species, atom, angular momentum
!  k, p, r, t for the 3-index tensor moment and vector component
integer, allocatable :: itmfix(:,:)
! 3-index tensor component with conventional normalisation
real(8), allocatable :: wkprfix(:)
! density matrices corresponding to the fixed tensor moments
complex(8), allocatable :: dmftm(:,:,:,:,:)
! fixed tensor moment potential matrix
complex(8), allocatable :: vmftm(:,:,:,:,:)
! fixed tensor moment step size
real(8) tauftm
! tensor moments at t=0 of a TDDFT+U calculations
real(8), allocatable :: wkpr0(:,:)
! tm3old is .true. if the tensor moments should be in the old complex convention
logical tm3old

end module

