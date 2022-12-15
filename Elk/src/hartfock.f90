
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hartfock
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,lp,nthd
real(8) etp,de
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:)
real(8), allocatable :: bmt(:,:,:),bir(:,:)
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
call init2
! read the charge density from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! generate the core wavefunctions and densities
call gencore
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW and local-orbital radial functions and integrals
call genapwlofr
! generate the spin-orbit coupling radial functions
call gensocfr
! generate the first- and second-variational eigenvectors and eigenvalues
call genevfsv
! find the occupation numbers and Fermi energy
call occupy
! generate the kinetic matrix elements in the first-variational basis
call genkmat(.true.,.true.)
! allocate local arrays
allocate(vmt(npcmtmax,natmtot),vir(ngtot))
if (hybrid.and.spinpol) then
  allocate(bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag))
end if
! only the MPI master process should write files
if (mp_mpi) then
! open HF_INFO.OUT file
  open(60,file='HF_INFO'//trim(filext),form='FORMATTED')
! open TOTENERGY.OUT
  open(61,file='TOTENERGY'//trim(filext),form='FORMATTED')
! open FERMIDOS.OUT
  open(62,file='FERMIDOS'//trim(filext),form='FORMATTED')
! open MOMENT.OUT if required
  if (spinpol) open(63,file='MOMENT'//trim(filext),form='FORMATTED')
! open GAP.OUT
  open(64,file='GAP'//trim(filext),form='FORMATTED')
! open DTOTENERGY.OUT
  open(66,file='DTOTENERGY'//trim(filext),form='FORMATTED')
! open MOMENTM.OUT
  if (spinpol) open(68,file='MOMENTM'//trim(filext),form='FORMATTED')
! write out general information to HF_INFO.OUT
  call writeinfo(60)
end if
! set last self-consistent loop flag
tlast=.false.
etp=0.d0
! begin the self-consistent loop
if (mp_mpi) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop started |")')
  write(60,'("+------------------------------+")')
end if
do iscl=1,maxscl
  if (mp_mpi) then
    write(60,*)
    write(60,'("+--------------------+")')
    write(60,'("| Loop number : ",I4," |")') iscl
    write(60,'("+--------------------+")')
    flush(60)
    write(*,*)
    write(*,'("Info(hartfock): self-consistent loop number : ",I4)') iscl
  end if
  if (iscl.ge.maxscl) then
    if (mp_mpi) then
      write(60,*)
      write(60,'("Reached self-consistent loops maximum")')
    end if
    tlast=.true.
  end if
! reset the OpenMP thread variables
  call omp_reset
! compute the Hartree-Fock local potentials
  call hflocal(vmt,vir,bmt,bir)
! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)
  call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsv) &
!$OMP NUM_THREADS(nthd)
  allocate(evecsv(nstsv,nstsv))
!$OMP DO
  do ik=1,nkpt
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
    call getevecsv(filext,ik,vkl(:,ik),evecsv)
! solve the Hartree-Fock eigenvalue equation
    call eveqnhf(ik,vmt,vir,bmt,bir,evecsv)
! write the eigenvalues/vectors to file
    call putevalsv(filext,ik,evalsv(:,ik))
    call putevecsv(filext,ik,evecsv)
  end do
!$OMP END DO
  deallocate(evecsv)
!$OMP END PARALLEL
  call freethd(nthd)
! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)
! broadcast eigenvalue array to every process
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(evalsv(:,ik),nstsv,mpi_double_precision,lp,mpicom,ierror)
  end do
! find the occupation numbers and Fermi energy
  call occupy
  if (mp_mpi) then
! write the occupation numbers to file
    do ik=1,nkpt
      call putoccsv(filext,ik,occsv(:,ik))
    end do
! write out the eigenvalues and occupation numbers
    call writeeval
! write the Fermi energy to file
    call writefermi
  end if
! generate the density and magnetisation
  call rhomag
! compute the energy components
  call energy
  if (mp_mpi) then
! output energy components
    call writeengy(60)
    write(60,*)
    write(60,'("Density of states at Fermi energy : ",G18.10)') fermidos
    write(60,'(" (states/Hartree/unit cell)")')
    write(60,*)
    write(60,'("Estimated indirect band gap : ",G18.10)') bandgap(1)
    write(60,'(" from k-point ",I6," to k-point ",I6)') ikgap(1),ikgap(2)
    write(60,'("Estimated direct band gap   : ",G18.10)') bandgap(2)
    write(60,'(" at k-point ",I6)') ikgap(3)
! write total energy to TOTENERGY.OUT
    write(61,'(G22.12)') engytot
    flush(61)
! write DOS at Fermi energy to FERMIDOS.OUT
    write(62,'(G18.10)') fermidos
    flush(62)
! output charges and moments
    call writechg(60)
    if (spinpol) then
      call writemom(60)
! write total moment to MOMENT.OUT
      write(63,'(3G18.10)') momtot(1:ndmag)
      flush(63)
! write total moment magnitude to MOMENTM.OUT
      write(68,'(G18.10)') momtotm
      flush(68)
    end if
! write estimated Hartree-Fock indirect band gap
    write(64,'(G22.12)') bandgap(1)
    flush(64)
  end if
  if (tlast) goto 10
! compute the change in total energy and check for convergence
  if (iscl.ge.2) then
    de=abs(engytot-etp)
    if (mp_mpi) then
      write(60,*)
      write(60,'("Absolute change in total energy (target) : ",G18.10," (",&
       &G18.10,")")') de,epsengy
    end if
    if (de.lt.epsengy) then
      if (mp_mpi) then
        write(60,*)
        write(60,'("Energy convergence target achieved")')
      end if
      tlast=.true.
    end if
    if (mp_mpi) then
      write(66,'(G18.10)') de
      flush(66)
    end if
  end if
  etp=engytot
! check for STOP file
  call checkstop
  if (tstop) tlast=.true.
! broadcast tlast from master process to all other processes
  call mpi_bcast(tlast,1,mpi_logical,0,mpicom,ierror)
end do
10 continue
if (mp_mpi) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop stopped |")')
  write(60,'("+------------------------------+")')
  if (maxscl.gt.1) then
    call writestate
    write(60,*)
    write(60,'("Wrote STATE.OUT")')
  end if
end if
!-----------------------!
!     compute forces    !
!-----------------------!
if (tforce) then
  call force
! output forces to HF_INFO.OUT
  if (mp_mpi) call writeforces(60)
end if
if (mp_mpi) then
  write(60,*)
  write(60,'("+----------------------------+")')
  write(60,'("| Elk version ",I1.1,".",I1.1,".",I2.2," stopped |")') version
  write(60,'("+----------------------------+")')
! close the HF_INFO.OUT file
  close(60)
! close the TOTENERGY.OUT file
  close(61)
! close the FERMIDOS.OUT file
  close(62)
! close the MOMENT.OUT and MOMENTM.OUT files
  if (spinpol) then
    close(63); close(68)
  end if
! close the GAP.OUT file
  close(64)
! close the DTOTENERGY.OUT file
  close(66)
end if
deallocate(vmt,vir)
if (hybrid.and.spinpol) deallocate(bmt,bir)
end subroutine

