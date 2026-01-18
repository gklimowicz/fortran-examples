
! Copyright (C) 2018 A. Davydov, P. Elliott, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gwscrho
use modmain
use modgw
use modmpi
use modomp
implicit none
! local variables
integer ik,nthd
integer nmix,nwork
real(8) dv
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:)
real(8), allocatable :: bmt(:,:,:),bir(:,:)
real(8), allocatable :: work(:)
complex(8), allocatable :: se(:,:,:)
! initialise universal variables
call init0
call init1
call init2
call init3
! read density and potentials from file
call readstate
call genvsig
! read the Fermi energy
call readfermi
! zero the exchange-correlation potential and magnetic field
vxcmt(:,:)=0.d0
vxcir(:)=0.d0
if (spinpol) then
  bxcmt(:,:,:)=0.d0
  bxcir(:,:)=0.d0
end if
! re-calculate the ground-state
call gencore
call linengy
call genapwlofr
call gensocfr
call genevfsv
call occupy
allocate(vmt(npcmtmax,natmtot),vir(ngtot))
if (spinpol) then
  allocate(bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag))
end if
if (mp_mpi) then
! open GW_INFO.OUT file
  open(60,file='GW_INFO.OUT',form='FORMATTED')
! open FERMIDOS.OUT
  open(62,file='FERMIDOS.OUT',form='FORMATTED')
! open MOMENT.OUT if required
  if (spinpol) open(63,file='MOMENT.OUT',form='FORMATTED')
! open GAP.OUT
  open(64,file='GAP.OUT',form='FORMATTED')
! open RMSDVS.OUT
  open(65,file='RMSDVS.OUT',form='FORMATTED')
! open MOMENTM.OUT
  if (spinpol) open(68,file='MOMENTM.OUT',form='FORMATTED')
! open RESIDUAL.OUT
  open(69,file='RESIDUAL.OUT',form='FORMATTED')
! write out general information to GW_INFO.OUT
  call writeigw(60)
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop started |")')
  write(60,'("+------------------------------+")')
end if
! size of mixing vector
nmix=size(vsbs)
! determine the size of the mixer work array
nwork=-1
call mixerifc(mixtype,nmix,vsbs,dv,nwork,vsbs)
allocate(work(nwork))
! initialise the mixer
iscl=0
call mixerifc(mixtype,nmix,vsbs,dv,nwork,work)
! set the stop signal to .false.
tstop=.false.
! set last self-consistent loop flag
tlast=.false.
do iscl=1,maxscl
  if (mp_mpi) then
    write(60,*)
    write(60,'("+--------------------+")')
    write(60,'("| Loop number : ",I4," |")') iscl
    write(60,'("+--------------------+")')
  end if
  call genpmat
  call epsinv
! compute the matrix elements of -V_xc and -B_xc
  call gwlocal(vmt,vir,bmt,bir)
! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)
  if (mp_mpi) write(*,*)
! loop over reduced k-point set
  call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(se) &
!$OMP NUM_THREADS(nthd)
  allocate(se(nstsv,nstsv,0:nwfm))
!$OMP DO
  do ik=1,nkpt
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(gwscrho_)
    write(*,'("Info(gwscrho): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(gwscrho_)
! determine the self-energy at the fermionic frequencies for current k-point
    call gwsefmk(ik,vmt,vir,bmt,bir,se)
! write the self-energy to file
    call putgwsefm(ik,se)
  end do
!$OMP END DO
  deallocate(se)
!$OMP END PARALLEL
  call freethd(nthd)
! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)
  if (mp_mpi) then
    write(60,*)
    write(60,'("Kohn-Sham Fermi energy : ",G18.10)') efermi
  end if
! determine the GW Fermi energy
  call gwefermi
  if (mp_mpi) then
    write(60,'("GW Fermi energy        : ",G18.10)') efermi
    flush(60)
  end if
! determine the density and magnetisation
  call gwrhomag
! invert the Kohn-Sham equations to find V_s and B_s
  call ksinvert
! mix the old effective potential and field with the new
  call mixerifc(mixtype,nmix,vsbs,dv,nwork,work)
! Fourier transform Kohn-Sham potential to G-space
  call genvsig
  if (mp_mpi) then
! write the Kohn-Sham occupation numbers to file
    do ik=1,nkpt
      call putoccsv(filext,ik,occsv(:,ik))
    end do
    call writeeval
    call writefermi
! write STATE.OUT file
    call writestate
    write(60,*)
    write(60,'("Density of states at Fermi energy : ",G18.10)') fermidos
    write(60,'(" (states/Hartree/unit cell)")')
    write(60,*)
    write(60,'("Estimated indirect band gap : ",G18.10)') bandgap(1)
    write(60,'(" from k-point ",I6," to k-point ",I6)') ikgap(1),ikgap(2)
    write(60,'("Estimated direct band gap   : ",G18.10)') bandgap(2)
    write(60,'(" at k-point ",I6)') ikgap(3)
! output charges and moments
    call writechg(60)
    if (spinpol) call writemom(60)
    write(60,*)
    write(60,'("Magnitude of residual : ",G18.10)') resoep
    flush(60)
! write DOS at Fermi energy to FERMIDOS.OUT
    write(62,'(G18.10)') fermidos
    flush(62)
    if (spinpol) then
! write total moment to MOMENT.OUT
      write(63,'(3G18.10)') momtot(1:ndmag)
      flush(63)
! write total moment magnitude to MOMENTM.OUT
      write(68,'(G18.10)') momtotm
      flush(68)
    end if
! write estimated Kohn-Sham indirect band gap
    write(64,'(G22.12)') bandgap(1)
    flush(64)
! write residual to RESIDUAL.OUT
    write(69,'(G18.10)') resoep
    flush(69)
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
! end the self-consistent loop
end do
10 continue
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
if (mp_mpi) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop stopped |")')
  write(60,'("+------------------------------+")')
! close the GW_INFO.OUT file
  close(60)
! close the FERMIDOS.OUT file
  close(62)
! close the MOMENT.OUT and MOMENTM.OUT files
  if (spinpol) then
    close(63); close(68)
  end if
! close the GAP.OUT file
  close(64)
! close the RMSDVS.OUT file
  close(65)
! close the RESIDUAL.OUT file
  close(69)
end if
deallocate(vmt,vir)
if (spinpol) deallocate(bmt,bir)
end subroutine

