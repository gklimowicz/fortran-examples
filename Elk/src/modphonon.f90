
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modphonon
use modmain

!--------------------------!
!     phonon variables     !
!--------------------------!
! number of phonon branches (3*natmtot)
integer nbph
! current phonon q-point, species, atom and polarisation index
integer iqph,isph,iaph,iasph,ipph
! tphq0 is .true. if q = 0
logical tphq0
! number of vectors for writing out frequencies and eigenvectors
integer nphwrt
! vectors in lattice coordinates for writing out frequencies and eigenvectors
real(8), allocatable :: vqlwrt(:,:)
! Coulomb pseudopotential
real(8) mustar
! number of temperatures for the Eliashberg equations and thermal properties
integer ntemp
! phonon frequencies for all q-points
real(8), allocatable :: wphq(:,:)

!-----------------------------!
!     supercell variables     !
!-----------------------------!
! number of primitive unit cells in phonon supercell
integer nscph
! Cartesian offset vectors for each primitive cell in the supercell
real(8), allocatable :: vscph(:,:)
! phonon displacement distance
real(8) deltaph

!---------------------!
!     k+q-vectors     !
!---------------------!
! k+q-vectors in lattice coordinates
real(8), allocatable :: vkql(:,:)
! k+q-vectors in Cartesian coordinates
real(8), allocatable :: vkqc(:,:)

!------------------------------!
!     G+q-vector variables     !
!------------------------------!
! G+q-vectors in Cartesian coordinates
real(8), allocatable :: vgqc(:,:)
! G+q-vector lengths
real(8), allocatable :: gqc(:)
! regularised Coulomb Green's function in G+q-space
real(8), allocatable :: gclgq(:)
! spherical Bessel functions j_l(|G+q|R_mt)
real(8), allocatable :: jlgqrmt(:,:,:)
! spherical harmonics for G+q-vectors
complex(8), allocatable :: ylmgq(:,:)
! structure factors for G+q-vectors
complex(8), allocatable :: sfacgq(:,:)
! smooth step function form factors for all species and G+q-vectors
real(8), allocatable :: ffacgq(:,:)
! characteristic function derivative in G- and G+q-space
complex(8), allocatable :: dcfunig(:)
! characteristic function derivative in real-space
complex(8), allocatable :: dcfunir(:)

!--------------------------------!
!     G+k+q-vector variables     !
!--------------------------------!
! number of G+k+q-vector for each k-point
integer, allocatable :: ngkq(:,:)
! index from G+k+q-vectors to G-vectors
integer, allocatable :: igkqig(:,:,:)
! G+k+q-vectors in lattice and Cartesian coordinates
real(8), allocatable :: vgkql(:,:,:,:),vgkqc(:,:,:,:)
! G+k+q-vector lengths
real(8), allocatable :: gkqc(:,:,:)
! structure factors for the G+k+q-vectors
complex(8), allocatable :: sfacgkq(:,:,:,:)

!----------------------------------------------------------!
!     density functional perturbation theory variables     !
!----------------------------------------------------------!
! density derivative
complex(8), allocatable :: drhomt(:,:),drhoir(:)
! magnetisation derivative
complex(8), allocatable :: dmagmt(:,:,:),dmagir(:,:)
! Coulomb potential derivative
complex(8), allocatable :: dvclmt(:,:),dvclir(:)
! if tphdyn is .true. then the phonon dynamical matrix is being calculated
logical :: tphdyn=.false.
! nuclear potential without the self-term; used for the phonon dynamical matrix
complex(8), allocatable :: zvnmt(:)
! gradient of vsmt for the displaced muffin-tin
complex(8), allocatable :: gvsmt(:)
! combined target array for dvsmt, dvsir, dbsmt and dbsir
complex(8), allocatable, target :: dvsbs(:)
! Kohn-Sham potential derivative
complex(8), pointer :: dvsmt(:,:),dvsir(:)
! Kohn-Sham effective magnetic field derivative
complex(8), pointer :: dbsmt(:,:,:),dbsir(:,:)
! G+q-space interstitial Kohn-Sham potential derivative
complex(8), allocatable :: dvsig(:)
! spin-orbit coupling radial function derivative
complex(8), allocatable :: dsocfr(:,:)
! APW-APW Hamiltonian integral derivatives
complex(8), allocatable :: dhaa(:,:,:,:,:,:)
! local-orbital-APW Hamiltonian integral derivatives
complex(8), allocatable :: dhloa(:,:,:,:,:)
! local-orbital-local-orbital Hamiltonian integral derivatives
complex(8), allocatable :: dhlolo(:,:,:,:)
! real Gaunt coefficient array
real(8), allocatable :: gntyyy(:,:,:)
! smallest allowed perturbation theory denominator for eigenvector derivatives
real(8) epsdev
! Fermi energy derivative
real(8) defermi
! first-variational eigenvalue derivatives
real(8), allocatable :: devalfv(:,:,:)
! second-variational eigenvalue derivatives
real(8), allocatable :: devalsv(:,:)
! second-variational occupation number derivatives
real(8), allocatable :: doccsv(:,:)

!-------------------------------------------+
!     electron-phonon mean-field theory     |
!-------------------------------------------+
! energy change of the electron system
real(8) dengye
! energy change of the phonon system
real(8) dengyph
! sum of energy changes in electron and phonon systems
real(8) dengy
! phonon frequency cut-off below which modes are neglected
real(8) wphcut
! scale factor of the electron-phonon term and its mixing parameter
real(8) ephscf(2)
! anomalous is .true. if only the anomalous density matrix is to be used in the
! construction of the electron-phonon Hamiltonian
logical anomalous
! tephde is .true. if D = D0 + E, otherwise D = D0
logical tephde
! single-precision electron-phonon matrix element array for each k- and q-point
complex(4), allocatable :: ephmkq(:,:,:,:,:)

end module

