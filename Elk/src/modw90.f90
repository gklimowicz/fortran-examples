
! Copyright (C) 2015 Jon Lafuente and Manh Duc Le; 2017-18 Arsenii Gerasimov,
! Yaroslav Kvashnin and Lars Nordstrom. This file is distributed under the terms
! of the GNU General Public License. See the file COPYING for license details.

module modw90

!---------------------------------------!
!     Wannier90 interface variables     !
!---------------------------------------!
! seedname for all Wannier90 files
character(256) seedname
! number of extra lines to write to .win file
integer nxlwin
! extra lines to write to .win file
character(256), allocatable :: xlwin(:)
! number of Wannier functions to calculate
integer num_wann
! number of bands to pass to Wannier90
integer num_bands
! index to bands
integer, allocatable :: idxw90(:)
! number of iterations for the minimisation of omega
integer num_iter
! trial step for the line search minimisation
real(8) trial_step
! maximum number of nearest neighbours per k-point
integer, parameter :: num_nnmax=12
! total number of nearest neighbours for each k-point
integer nntot
! list of nearest neighbours for each k-point
integer, allocatable :: nnlist(:,:)
! G-vector offset for each nearest neighbour
integer, allocatable :: nncell(:,:,:)

end module

