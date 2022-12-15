
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine fermisurf
use modmain
use modomp
implicit none
! local variables
integer ik,nst,ist
integer ist0,ist1,jst0,jst1
integer i1,i2,i3,nf,f
integer np,i,nthd
real(8) e0,e1,prd,v(3)
! allocatable arrays
integer, allocatable :: idx(:)
real(8), allocatable :: evalfv(:,:),e(:),vpc(:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
! read density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW and local-orbital radial functions and integrals
call genapwlofr
! generate the spin-orbit coupling radial functions
call gensocfr
! begin parallel loop over reduced k-points set
call holdthd(nkpt,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv) &
!$OMP NUM_THREADS(nthd)
allocate(evalfv(nstfv,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
!$OMP DO
do ik=1,nkpt
!$OMP CRITICAL(fermisurf_)
  write(*,'("Info(fermisurf): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(fermisurf_)
! solve the first- and second-variational eigenvalue equations
  call eveqn(ik,evalfv,evecfv,evecsv)
! end loop over reduced k-points set
end do
!$OMP END DO
deallocate(evalfv,evecfv,evecsv)
!$OMP END PARALLEL
call freethd(nthd)
! if iterative diagonalisation is used the eigenvalues must be reordered
if (tefvit.and.(.not.spinpol)) then
  allocate(idx(nstsv),e(nstsv))
  do ik=1,nkpt
    e(:)=evalsv(:,ik)
    call sortidx(nstsv,e,idx)
    evalsv(1:nstsv,ik)=e(idx(1:nstsv))
  end do
  deallocate(idx,e)
end if
! generate the plotting point grid in Cartesian coordinates (this has the same
! arrangement as the k-point grid)
np=np3d(1)*np3d(2)*np3d(3)
allocate(vpc(3,np))
call plotpt3d(vpc)
do i=1,np
  v(:)=vpc(:,i)
  call r3mv(bvec,v,vpc(:,i))
end do
! number of files to plot (2 for collinear magnetism, 1 otherwise)
if (ndmag.eq.1) then
  nf=2
else
  nf=1
end if
do f=1,nf
  if (nf.eq.2) then
    if (f.eq.1) then
      open(50,file='FERMISURF_UP.OUT',form='FORMATTED',action='WRITE')
      jst0=1; jst1=nstfv
    else
      open(50,file='FERMISURF_DN.OUT',form='FORMATTED',action='WRITE')
      jst0=nstfv+1; jst1=2*nstfv
    end if
  else
    open(50,file='FERMISURF.OUT',form='FORMATTED',action='WRITE')
    jst0=1; jst1=nstsv
  end if
! find the range of eigenvalues which contribute to the Fermi surface (Lars)
  ist0=jst1; ist1=jst0
  do ist=jst0,jst1
    e0=minval(evalsv(ist,:)); e1=maxval(evalsv(ist,:))
! determine if the band crosses the Fermi energy
    if ((e0.lt.efermi).and.(e1.gt.efermi)) then
      ist0=min(ist0,ist); ist1=max(ist1,ist)
    end if
  end do
  nst=ist1-ist0+1
  if (task.eq.100) then
! write product of eigenstates minus the Fermi energy
    write(50,'(3I6," : grid size")') np3d(:)
    i=0
    do i3=0,ngridk(3)-1
      do i2=0,ngridk(2)-1
        do i1=0,ngridk(1)-1
          i=i+1
          ik=ivkik(i1,i2,i3)
          prd=product(evalsv(ist0:ist1,ik)-efermi)
          write(50,'(4G18.10)') vpc(:,i),prd
        end do
      end do
    end do
  else
! write the eigenvalues minus the Fermi energy separately
    write(50,'(4I6," : grid size, number of states")') np3d(:),nst
    i=0
    do i3=0,ngridk(3)-1
      do i2=0,ngridk(2)-1
        do i1=0,ngridk(1)-1
          i=i+1
          ik=ivkik(i1,i2,i3)
          write(50,'(3G18.10,40F14.8)') vpc(:,i),evalsv(ist0:ist1,ik)-efermi
        end do
      end do
    end do
  end if
  close(50)
end do
write(*,*)
write(*,'("Info(fermisurf):")')
if (ndmag.eq.1) then
  write(*,'(" 3D Fermi surface data written to FERMISURF_UP.OUT and &
   &FERMISURF_DN.OUT")')
else
  write(*,'(" 3D Fermi surface data written to FERMISURF.OUT")')
end if
if (task.eq.100) then
  write(*,'(" in terms of the product of eigenvalues minus the Fermi energy")')
else
  write(*,'(" in terms of separate eigenvalues minus the Fermi energy")')
end if
deallocate(vpc)
end subroutine

