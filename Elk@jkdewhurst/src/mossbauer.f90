
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: mossbauer
! !INTERFACE:
subroutine mossbauer
! !USES:
use modmain
use modmpi
use modtest
! !DESCRIPTION:
!   Computes the contact charge density and magnetic hyperfine field for each
!   atom and outputs the data to the file {\tt MOSSBAUER.OUT}.
!   See S. Bl\"{u}gel, H. Akai, R. Zeller, and P. H. Dederichs,
!   {\it Phys. Rev. B} {\bf 35}, 3271 (1987).
!
! !REVISION HISTORY:
!   Created May 2004 (JKD)
!   Contact hyperfine field evaluated at the nuclear radius rather than averaged
!   over the Thomson sphere, June 2019 (JKD)
!   Added spin and orbital dipole terms, July 2019 (JKD)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ia,ias
integer nr,nri,nrn
real(8) mc(3),bc(3),bd(3)
real(8) rho0,rhon,rhoa
real(8) cb,t1
! allocatable arrays
real(8), allocatable :: fr(:)
! spin dipole field prefactor
cb=gfacte/(4.d0*solsc)
! initialise universal variables
call init0
call init1
! read density and potentials from file
call readstate
! generate the core wavefunctions and densities
call gencore
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW and local-orbital radial functions and integrals
call genapwlofr
! read the eigenvalues and occupancies from file
call readevalsv
call readoccsv
! calculate the density
call rhomag
! calculate the dipole magnetic field if required
if (tbdip) call bdipole
! allocate local arrays
allocate(fr(nrmtmax))
if (mp_mpi) then
  open(50,file='MOSSBAUER.OUT',form='FORMATTED')
end if
! loop over species
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  nrn=nrnucl(is)
! loop over atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    if (mp_mpi) then
      write(50,*)
      write(50,*)
      write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),&
       ia
      write(50,*)
      write(50,'(" approximate nuclear radius : ",G18.10)') rnucl(is)
      write(50,'(" number of mesh points to nuclear radius : ",I6)') nrn
    end if
!--------------------------------!
!     contact charge density     !
!--------------------------------!
! extract the l=m=0 component of the muffin-tin density
    call rfmtlm(1,nr,nri,rhomt(:,ias),fr)
    rho0=fr(1)*y00
    rhon=fr(nrn)*y00
    t1=dot_product(wrmt(1:nrn,is),fr(1:nrn))
    rhoa=fourpi*y00*t1/volnucl(is)
    if (mp_mpi) then
      write(50,*)
      write(50,'(" density at nuclear center      : ",G18.10)') rho0
      write(50,'(" density at nuclear surface     : ",G18.10)') rhon
      write(50,'(" average contact charge density : ",G18.10)') rhoa
    end if
!----------------------------------!
!     magnetic hyperfine field     !
!----------------------------------!
    if (spinpol) then
! contact term
      do idm=1,ndmag
! extract the l=m=0 component of the muffin-tin magnetisation
        call rfmtlm(1,nr,nri,magmt(:,ias,idm),fr)
        t1=dot_product(wrmt(1:nrn,is),fr(1:nrn))
        mc(idm)=fourpi*y00*t1/volnucl(is)
      end do
      t1=8.d0*pi*cb/(3.d0*solsc)
      bc(1:ndmag)=t1*mc(1:ndmag)
      if (mp_mpi) then
        write(50,*)
        write(50,'(" contact magnetic moment (mu_B) : ",3G18.10)') mc(1:ndmag)
        write(50,'(" contact field : ",3G18.10)') bc(1:ndmag)
        write(50,'("  tesla : ",3G18.10)') b_si*bc(1:ndmag)
      end if
! spin and orbital dipole term
      if (tbdip) then
! extract the l=m=0 component of the dipole field
        do idm=1,3
          call rfmtlm(1,nr,nri,bdmt(:,ias,idm),fr)
          t1=dot_product(wrmt(1:nrn,is),fr(1:nrn))
          bd(idm)=fourpi*y00*t1/(volnucl(is)*solsc)
        end do
        if (mp_mpi) then
          write(50,*)
          if (tjr) then
            write(50,'(" spin and orbital dipole field : ",3G18.10)') bd
          else
            write(50,'(" spin dipole field : ",3G18.10)') bd
          end if
          write(50,'("  tesla : ",3G18.10)') b_si*bd
        end if
! write to test file if required
        call writetest(110,'hyperfine field',nv=3,tol=1.d-4,rva=bd)
      end if
    end if
  end do
end do
if (mp_mpi) then
  if (spinpol.and.tbdip) then
    write(50,*)
    write(50,'("Note that the contact term is implicitly included in the &
     &spin dipole field")')
    write(50,'(" but may not match exactly with the directly &
     &calculated value.")')
  end if
  close(50)
  write(*,*)
  write(*,'("Info(mossbauer):")')
  write(*,'(" Mossbauer parameters written to MOSSBAUER.OUT")')
end if
deallocate(fr)
end subroutine
!EOC

