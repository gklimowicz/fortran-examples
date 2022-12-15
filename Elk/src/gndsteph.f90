
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gndsteph
use modmain
use modphonon
use modbog
use modmpi
use modomp
implicit none
! local variables
integer nmix,nwork
real(8) dv
! allocatable arrays
real(8), allocatable :: work(:)
! initialise universal variables
call init0
call init1
call init2
call readstate
call genvsig
call gencore
call readfermi
call linengy
call genapwlofr
call gensocfr
call genevfsv
! precise determination of the Fermi energy
swidth0=swidth
swidth=1.d-5
call occupy
swidth=swidth0
! initialise electron-phonon variables
call initeph
! size of mixing vector for electron and phonon density matrices (complex array)
nmix=2*size(duvwx)
! determine the size of the mixer work array
nwork=-1
call mixerifc(mixtype,nmix,duvwx,dv,nwork,duvwx)
allocate(work(nwork))
! initialise the mixer
iscl=0
call mixerifc(mixtype,nmix,duvwx,dv,nwork,work)
! set the stop signal to .false.
tstop=.false.
! set last self-consistent loop flag
tlast=.false.
! only the MPI master process should write files
if (mp_mpi) then
! open EPH_INFO.OUT file
  open(60,file='EPH_INFO.OUT',form='FORMATTED')
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop started |")')
  write(60,'("+------------------------------+")')
! open EPHGAP.OUT
  open(64,file='EPHGAP.OUT',form='FORMATTED')
! open RMSDVS.OUT
  open(65,file='RMSDVS.OUT',form='FORMATTED')
! open FACE.OUT
  open(67,file='FACE.OUT',form='FORMATTED')
end if
if (mp_mpi) write(*,*)
! begin the self-consistent loop
do iscl=1,maxscl
  if (mp_mpi) then
    write(60,*)
    write(60,'("+--------------------+")')
    write(60,'("| Loop number : ",I4," |")') iscl
    write(60,'("+--------------------+")')
    flush(60)
    write(*,'("Info(gndsteph): self-consistent loop number : ",I4)') iscl
  end if
  if (iscl.ge.maxscl) then
    if (mp_mpi) then
      write(60,*)
      write(60,'("Reached self-consistent loops maximum")')
    end if
    tlast=.true.
  end if
! determine change in electron and phonon energies
  call dengyeph
! solve the electron and phonon eigenvalue equations
  call eveqneph
! update the Fermi energy
  call occupyuv
  if (mp_mpi) then
! write the electronic eigenvalues to file
    call writeevaluv
! write the phononic eigenvalues to file
    call writeevalwx
! write the Fermi energy to file
    call writefermi
    write(60,*)
    write(60,'("Energies :")')
    write(60,'(" Fermi",T30,": ",G22.12)') efermi
    write(60,'(" electronic change",T30,": ",G22.12)') dengye
    write(60,'(" phononic change",T30,": ",G22.12)') dengyph
    write(60,'(" sum of changes",T30,": ",G22.12)') dengy
    write(60,*)
    write(60,'("Estimated indirect band gap : ",G18.10)') bandgap(1)
    write(60,'(" from k-point ",I6," to k-point ",I6)') ikgap(1),ikgap(2)
    write(60,'("Estimated direct band gap   : ",G18.10)') bandgap(2)
    write(60,'(" at k-point ",I6)') ikgap(3)
    write(60,*)
    write(60,'("Fermionic anomalous correlation entropy : ",G18.10)') face
    write(60,*)
    write(60,'("Electron-phonon scaling factor : ",G18.10)') ephscf(1)
! write estimated indirect band gap
    write(64,'(G22.12)') bandgap(1)
    flush(64)
! write the fermionic anomalous correlation entropy
    write(67,'(G18.10)') face
    flush(67)
  end if
! mix the old and new electron and phonon density matrices
  call mixerifc(mixtype,nmix,duvwx,dv,nwork,work)
! adjust the electron-phonon term scale factor towards 1
  ephscf(1)=(1.d0-ephscf(2))*ephscf(1)+ephscf(2)
! exit self-consistent loop if required
  if (tlast) goto 10
! check for convergence
  if (iscl.ge.2) then
    if (mp_mpi) then
      write(60,*)
      write(60,'("RMS change in density matrices (target) : ",G18.10," (",&
       &G18.10,")")') dv,epspot
      write(65,'(G18.10)') dv
      flush(65)
      if (dv.lt.epspot) then
        write(60,*)
        write(60,'("Convergence targets achieved")')
        tlast=.true.
      end if
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
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop stopped |")')
  write(60,'("+------------------------------+")')
! close the EPH_INFO.OUT file
  close(60)
! close the EPHGAP.OUT file
  close(64)
! close the RMSDVS.OUT file
  close(65)
! close the FACE.OUT file
  close(67)
end if
deallocate(work)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

