
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gndstulr
use modmain
use modulr
use moddftu
use modmpi
use modomp
implicit none
! local variables
logical twrite
integer ik0,ir,lp,nthd
integer nmix,nwork,n
real(8) dv
! allocatable arrays
integer(8), allocatable :: lock(:)
real(8), allocatable :: work(:)
complex(8), allocatable :: evecu(:,:)
if (xctype(1).lt.0) then
  write(*,*)
  write(*,'("Error(gndstulr): ultra long-range does not work with OEP")')
  write(*,*)
  stop
end if
if (spincore) then
  write(*,*)
  write(*,'("Error(gndstulr): ultra long-range does not work with &
   &spin-polarised cores")')
  write(*,*)
  stop
end if
! no k-point reduction
reducek0=reducek
reducek=0
! initialise global variables
call init0
call init1
! write the kappa-points to file
call writekpa
! write the k+kappa-points to file
call writekpts
! read the regular Kohn-Sham potential from file
call readstate
! generate the first- and second-variational eigenvectors and eigenvalues for
! the k+kappa-point set
call genvsig
call gencore
call readfermi
call linengy
call genapwlofr
call gensocfr
call genevfsv
call occupy
! initialise the ultra long-range variables
call initulr
if (task.eq.700) then
! initialise the long-range Kohn-Sham potential and magnetic field
  call potuinit
else
! read in the potential and density from STATE_ULR.OUT
  call readstulr
end if
! initialise the external Coulomb potential
call vclqinit
! write the long-range Hamiltonian diagonal blocks to file
call genhdbulr
! size of mixing vector (complex array)
nmix=2*size(vsbsq)
! determine the size of the mixer work array
nwork=-1
call mixerifc(mixtype,nmix,vsbsq,dv,nwork,vsbsq)
allocate(work(nwork))
! initialise the mixer
iscl=0
call mixerifc(mixtype,nmix,vsbsq,dv,nwork,work)
! initialise the OpenMP locks
allocate(lock(nqpt))
do ir=1,nqpt
  call omp_init_lock(lock(ir))
end do
! set last self-consistent loop flag
tlast=.false.
! begin the self-consistent loop
if (mp_mpi) then
! open ULR_INFO.OUT file
  open(60,file='ULR_INFO.OUT',form='FORMATTED')
! open RMSDVS.OUT
  open(65,file='RMSDVS.OUT',form='FORMATTED')
  call writeinfou(60)
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
  end if
  if (iscl.ge.maxscl) then
    if (mp_mpi) then
      write(60,*)
      write(60,'("Reached self-consistent loops maximum")')
    end if
    write(*,*)
    write(*,'("Warning(gndstulr): failed to reach self-consistency after ",I4,&
     &" loops")') iscl
    tlast=.true.
  end if
! reset the OpenMP thread variables
  call omp_reset
! apply required local operations to the potential and magnetic field
  call vblocalu
! zero the density and magnetisation
  rhormt(:,:,:)=0.d0
  rhorir(:,:)=0.d0
  if (spinpol) then
    magrmt(:,:,:,:)=0.d0
    magrir(:,:,:)=0.d0
  end if
! loop over original k-points
  call holdthd(nkpt0/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecu) &
!$OMP NUM_THREADS(nthd)
  allocate(evecu(nstulr,nstulr))
!$OMP DO
  do ik0=1,nkpt0
! distribute among MPI processes
    if (mod(ik0-1,np_mpi).ne.lp_mpi) cycle
! solve the ultra long-range eigenvalue equation
    call eveqnulr(ik0,evecu)
! add to the density, magnetisation and current
    call rhomaguk(ik0,lock,evecu)
  end do
!$OMP END DO
  deallocate(evecu)
!$OMP END PARALLEL
  call freethd(nthd)
  if (np_mpi.gt.1) then
! broadcast eigenvalue array to every process
    do ik0=1,nkpt0
      lp=mod(ik0-1,np_mpi)
      call mpi_bcast(evalu(:,ik0),nstulr,mpi_double_precision,lp,mpicom,ierror)
    end do
! add densities from each process and redistribute
    n=npcmtmax*natmtot*nqpt
    call mpi_allreduce(mpi_in_place,rhormt,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
    n=ngtc*nqpt
    call mpi_allreduce(mpi_in_place,rhorir,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
    if (spinpol) then
      n=npcmtmax*natmtot*ndmag*nqpt
      call mpi_allreduce(mpi_in_place,magrmt,n,mpi_double_precision,mpi_sum, &
       mpicom,ierror)
      n=ngtc*ndmag*nqpt
      call mpi_allreduce(mpi_in_place,magrir,n,mpi_double_precision,mpi_sum, &
       mpicom,ierror)
    end if
  end if
! find the occupation numbers and Fermi energy
  call occupyulr
! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)
! add the core density
  call rhocoreu
! perform partial Fourier transform to Q-space
  call rhomagq
! determine the muffin-tin and interstitial charges and moments
  call chargeu
  call momentu
! compute the ultra long-range Kohn-Sham potential
  call potksu
! mix the old potential and field with the new
  call mixerifc(mixtype,nmix,vsbsq,dv,nwork,work)
! multiply the RMS change in potential by the number of Q-points
  dv=dv*dble(nfqrz)
! calculate and add the fixed spin moment effective field (after mixing)
  call bfieldfsm
  call addbfsmu
  if (mp_mpi) then
! write eigenvalues to file
    call writeevalu
! output energy components
    call writeengyu(60)
! output charges
    call writechg(60)
! write muffin-tin charges for each R-vector
    call writechgrmt
    if (spinpol) then
! output moments
      call writemom(60)
! write muffin-tin moments for each R-vector
      call writemomrmt
    end if
! output effective fields for fixed spin moment calculations
    if (fsmtype.ne.0) call writefsm(60)
! check for existence of the WRITE file
    call checkwrite(twrite)
! check self-consistent loop number modulo nwrite
    if (nwrite.ge.1) then
      if (mod(iscl,nwrite).eq.0) twrite=.true.
    end if
! write STATE_ULR.OUT file if required
    if (twrite) then
      call writestulr
      write(60,*)
      write(60,'("Wrote STATE_ULR.OUT")')
    end if
  end if
! exit self-consistent loop if required
  if (tlast) goto 10
! check for convergence
  if (iscl.ge.2) then
    if (mp_mpi) then
      write(60,*)
      write(60,'("RMS change in Kohn-Sham potential (target) : ",G18.10," (",&
       &G18.10,")")') dv,epspot
      flush(60)
      write(65,'(G18.10)') dv
      flush(65)
    end if
    if (dv.lt.epspot) then
      if (mp_mpi) then
        write(60,*)
        write(60,'("Convergence targets achieved")')
      end if
      tlast=.true.
    end if
  end if
! check for STOP file
  call checkstop
  if (tstop) tlast=.true.
! broadcast tlast from master process to all other processes
  call mpi_bcast(tlast,1,mpi_logical,0,mpicom,ierror)
! reset the OpenMP thread variables
  call omp_reset
end do
10 continue
if (mp_mpi) then
! output timing information
  write(60,*)
  write(60,'("Timings (CPU seconds) :")')
  write(60,'(" Hamiltonian matrix set up",T40,": ",F12.2)') timemat
  write(60,'(" eigenvalue equation",T40,": ",F12.2)') timesv
  write(60,'(" charge density calculation",T40,": ",F12.2)') timerho
  write(60,'(" potential calculation",T40,": ",F12.2)') timepot
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop stopped |")')
  write(60,'("+------------------------------+")')
  if (maxscl.gt.1) then
    call writestulr
    write(60,*)
    write(60,'("Wrote STATE_ULR.OUT")')
  end if
! close the ULR_INFO.OUT file
  close(60)
! close the RMSDVS.OUT file
  close(65)
end if
! destroy the OpenMP locks
do ir=1,nqpt
  call omp_destroy_lock(lock(ir))
end do
deallocate(lock,work)
! restore original parameters
reducek=reducek0
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

