
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genapwlofr
use modomp
implicit none
! local variables
integer nthd
call holdthd(2,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTIONS
!$OMP SECTION
! generate the APW radial functions
call genapwfr
!$OMP SECTION
! generate the local-orbital radial functions
call genlofr
!$OMP END SECTIONS
!$OMP SECTIONS
!$OMP SECTION
! compute the overlap radial integrals
call olprad
!$OMP SECTION
! compute the Hamiltonian radial integrals
call hmlrad
!$OMP END SECTIONS
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

