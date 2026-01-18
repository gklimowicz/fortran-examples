
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modbog

! combined target array for fermionic and bosonic density matrices
complex(8), allocatable, target :: duvwx(:)

!----------------------------------------!
!     fermionic Bogoliubov variables     !
!----------------------------------------!
! Bogoliubov equation eigenvalues
real(8), allocatable :: evaluv(:,:)
! V-norm for each state and k-point
real(8), allocatable :: vnorm(:,:)
! Fermi energy adjustment step size
real(8) tauefm
! Fermi energy convergence tolerance
real(8) epsefm
! Hartree-Fock-Bogoliubov coupling constant
real(8) ehfb
! density matrices VV† and UV†
complex(8), pointer :: dvv(:,:,:),duv(:,:,:)
! fermionic anomalous correlation entropy
real(8) face
! bdiag is .true. if the matrix B is taken to be diagonal
logical bdiag
! cut-off energy for matrix B (elements outside this window are set to zero)
real(8) ecutb

!--------------------------------------!
!     bosonic Bogoliubov variables     !
!--------------------------------------!
! Bogoliubov equation eigenvalues
real(8), allocatable :: evalwx(:,:)
! power used in formula for (W,X) pseudo-normalisation (see article)
integer pwxpsn
! X-norm for each branch and q-point
real(8), allocatable :: xnorm(:,:)
! density matrices XX† and WX†
complex(8), pointer :: dxx(:,:,:),dwx(:,:,:)
! ediag is .true. if the matrix E is taken to be diagonal
logical ediag

end module

