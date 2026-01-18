
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writehmlbse
use modmain
use modmpi
! sets up the BSE matrix and writes it to file
implicit none
! local variables
integer ik,jk,ist,jst
integer a,b,i,j,m,n
real(8) t1
! initialise global variables
call init0
call init1
call init2
call init3
! read density and potentials from file
call readstate
! read Fermi energy from a file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! check if system is metallic
t1=minval(abs(0.5d0-occsv(:,:)/occmax))
if (abs(t1-0.5d0).gt.0.01d0) then
  write(*,*)
  write(*,'("Warning(writehmlbse): system is metallic, the BSE may fail")')
  write(*,'("Try using a different vkloff or reducing swidth")')
end if
! generate the BSE state index arrays
call genidxbse
if (allocated(hmlbse)) deallocate(hmlbse)
allocate(hmlbse(nmbse,nmbse))
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writehmlbse): setting up BSE Hamiltonian matrix")')
end if
! zero the BSE Hamiltonian
hmlbse(:,:)=0.d0
! compute diagonal matrix elements
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
  do i=1,nvbse
    ist=istbse(i,ik)
    do j=1,ncbse
      jst=jstbse(j,ik)
      a=ijkbse(i,j,ik)
      hmlbse(a,a)=(evalsv(jst,jk)+scissor)-evalsv(ist,jk)
      if (bsefull) then
        b=a+nbbse
        hmlbse(b,b)=-hmlbse(a,a)
      end if
    end do
  end do
end do
! add the exchange matrix elements
if (hxbse) call hmlxbse
! add the direct matrix elements
if (hdbse) call hmldbse
! add matrices from all processes and redistribute
if (np_mpi.gt.1) then
! ensure that the number of elements transmitted by MPI is not larger than the
! maximum packet size (assumed to be 1024³ bytes)
  m=67108864/nmbse
  do b=1,nmbse,m
    n=min(nmbse-b+1,m)
    n=nmbse*n
    call mpi_allreduce(mpi_in_place,hmlbse(1,b),n,mpi_double_complex,mpi_sum, &
     mpicom,ierror)
  end do
end if
! write the BSE matrix to HMLBSE.OUT
if (mp_mpi) then
  open(140,file='HMLBSE.OUT',form='UNFORMATTED',action='WRITE')
  write(140) nmbse
  write(140) hmlbse
  close(140)
  write(*,*)
  write(*,'("Info(writehmlbse): BSE Hamiltonian matrix written to HMLBSE.OUT")')
end if
! deallocate global BSE arrays
deallocate(istbse,jstbse,ijkbse,hmlbse)
end subroutine

