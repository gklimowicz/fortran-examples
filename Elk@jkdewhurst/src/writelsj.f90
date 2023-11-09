
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl and
! F. Cricchio. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine writelsj
use modmain
use moddftu
use modmpi
use modtest
implicit none
! local variables
integer kst,ik,ist
integer is,ia,ias
real(8) xl(3),xs(3)
! allocatable arrays
real(8), allocatable :: xj(:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
! initialise universal variables
call init0
call init1
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! allocate local arrays
allocate(xj(3,natmtot))
allocate(dmat(lmmaxo,nspinor,lmmaxo,nspinor,natmtot))
if (task.eq.15) then
! compute total L, S and J
! read in the occupation numbers
  do ik=1,nkpt
    call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
  end do
! generate the density matrix in each muffin-tin
  call gendmat(.false.,.false.,0,lmaxo,lmmaxo,dmat)
  if (mp_mpi) then
    open(50,file='LSJ.OUT',form='FORMATTED',action='WRITE')
    write(50,*)
    write(50,'("Expectation values are computed only over the muffin-tin")')
! loop over species and atoms
    do is=1,nspecies
      write(50,*)
      write(50,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! calculate the expectation value of L and S
        call dmatls(dmat(:,:,:,:,ias),xl,xs)
! J = L + S
        xj(:,ias)=xl(:)+xs(:)
        write(50,'(" atom : ",I4)') ia
        write(50,'("  L : ",3G18.10)') xl(:)
        write(50,'("  S : ",3G18.10)') xs(:)
        write(50,'("  J : ",3G18.10)') xj(:,ias)
! end loop over atoms and species
      end do
    end do
    close(50)
    write(*,*)
    write(*,'("Info(writelsj):")')
    write(*,'(" muffin-tin L, S and J expectation values written to LSJ.OUT")')
  end if
! write J to test file
  call writetest(15,'total muffin-tin angular momentum',nv=3*natmtot,tol=1.d-3,&
   rva=xj)
else
! compute L, S and J for all states in kstlist
  if (mp_mpi) then
    open(50,file='LSJ_KST.OUT',form='FORMATTED',action='WRITE')
    write(50,*)
    write(50,'("Expectation values are computed only over the muffin-tin")')
  end if
  do kst=1,nkstlist
    ik=kstlist(1,kst)
    ist=kstlist(2,kst)
    if ((ik.le.0).or.(ik.gt.nkpt)) then
      write(*,*)
      write(*,'("Error(writelsj): k-point out of range : ",I8)') ik
      write(*,*)
      stop
    end if
    if ((ist.le.0).or.(ist.gt.nstsv)) then
      write(*,*)
      write(*,'("Error(writelsj): state out of range : ",I8)') ist
      write(*,*)
      stop
    end if
! select a particular wavefunction using its occupancy
    occsv(:,:)=0.d0
    occsv(ist,ik)=1.d0/wkpt(ik)
! no symmetrisation required
    nsymcrys=1
! generate the density matrix in each muffin-tin
    call gendmat(.false.,.false.,0,lmaxo,lmmaxo,dmat)
    if (mp_mpi) then
! loop over species and atoms
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
! calculate the expectation value of L and S
          call dmatls(dmat(:,:,:,:,ias),xl,xs)
! J = L + S
          xj(:,ias)=xl(:)+xs(:)
          write(50,*)
          write(50,'("k-point : ",I6,3G18.10)') ik,vkl(:,ik)
          write(50,'("state : ",I6)') ist
          write(50,'("species : ",I4," (",A,"), atom : ",I4)') is, &
           trim(spsymb(is)),ia
          write(50,'(" L : ",3G18.10)') xl(:)
          write(50,'(" S : ",3G18.10)') xs(:)
          write(50,'(" J : ",3G18.10)') xj(:,ias)
        end do
      end do
    end if
  end do
  if (mp_mpi) then
    close(50)
    write(*,*)
    write(*,'("Info(writelsj):")')
    write(*,'(" muffin-tin L, S and J expectation values for each k-point and &
     &state")')
    write(*,'("  in kstlist written to LSJ_KST.OUT")')
  end if
! write J to test file
  call writetest(16,'muffin-tin angular momentum for one state',nv=3*natmtot, &
   tol=1.d-3,rva=xj)
end if
deallocate(xj,dmat)
end subroutine

