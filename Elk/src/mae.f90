
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mae
use modmain
use modmpi
use moddelf
use modtest
implicit none
! local variables
integer i,j,i0,i1
real(8) e0,e1,de
real(8) v1(3),v2(3),th
real(8) a(3,3),b(3,3)
! initialise global variables
call init0
! store original parameters
avec0(:,:)=avec(:,:)
spinpol0=spinpol
spinorb0=spinorb
cmagz0=cmagz
bfieldc00(:)=bfieldc0(:)
reducebf0=reducebf
fsmtype0=fsmtype
vkloff0(:)=vkloff(:)
! enable spin-orbit coupling
spinorb=.true.
! enforce collinear magnetism in the z-direction
cmagz=.true.
! no fixed spin moment calculation: the crystal is rotated instead
fsmtype=0
! if task=28 then start from atomic densities; if task=29 read STATE.OUT
if (task.eq.28) then
  trdstate=.false.
else
  trdstate=.true.
end if
! zero k-point offset
vkloff(:)=0.d0
! start with large magnetic field
bfieldc0(1:2)=0.d0
bfieldc0(3)=-2.d0
! reduce the external magnetic field after each s.c. loop
reducebf=0.75d0
! generate the spin moment directions in (theta,phi) coordinates
call gentpmae
! open MAE_INFO.OUT
if (mp_mpi) then
  open(71,file='MAE_INFO.OUT',form='FORMATTED')
  write(71,*)
  write(71,'("Scale factor of spin-orbit coupling term : ",G18.10)') socscf
end if
i0=1; i1=1
e0=1.d8; e1=-1.d8
! loop over points on sphere
do i=1,npmae
  if (mp_mpi) then
    write(*,'("Info(mae): fixed spin moment direction ",I6," of ",I6)') i,npmae
  end if
! rotate lattice vectors instead of moment (thanks to J. Glasbrenner,
! K. Bussmann and I. Mazin)
! first by -phi around the z-axis
  v1(:)=0.d0
  v1(3)=1.d0
  th=-tpmae(2,i)
  call axangrot(v1,th,a)
! then by -theta around the y-axis
  v1(:)=0.d0
  v1(2)=1.d0
  th=-tpmae(1,i)
  call axangrot(v1,th,b)
  call r3mm(b,a,rotsht)
  call r3mm(rotsht,avec0,avec)
! find the corresponding moment direction vector
  call r3minv(rotsht,a)
  v1(:)=0.d0
  v1(3)=1.d0
  call r3mv(a,v1,v2)
  do j=1,3
    if (abs(v2(j)).lt.epslat) v2(j)=0.d0
  end do
! rotate the spherical cover used for the spherical harmonic transform
  trotsht=.true.
! run the ground-state calculation
  call gndstate
! subsequent calculations should read the previous density
  trdstate=.true.
! make external magnetic field small
  bfieldc0(3)=-0.01d0
  if (mp_mpi) then
    write(71,*)
    write(71,'("Fixed spin moment direction point ",I6," of ",I6)') i,npmae
    write(71,'("Spherical coordinates of direction : ",2G18.10)') tpmae(:,i)
    write(71,'("Direction vector (Cartesian coordinates) : ",3G18.10)') v2
    write(71,'("Calculated total moment magnitude : ",G18.10)') momtotm
    write(71,'("Total energy : ",G22.12)') engytot
    flush(71)
  end if
! check for minimum and maximum total energy
  if (engytot.lt.e0) then
    e0=engytot
    i0=i
  end if
  if (engytot.gt.e1) then
    e1=engytot
    i1=i
  end if
! delete the eigenvector files
  call delfiles(evec=.true.)
! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)
end do
! magnetic anisotropy energy
de=e1-e0
if (mp_mpi) then
  write(71,*)
  write(71,'("Minimum energy point : ",I6)') i0
  write(71,'("Maximum energy point : ",I6)') i1
  write(71,*)
  write(71,'("Estimated magnetic anisotropy energy (MAE) : ",G18.10)') de
  write(71,*)
  write(71,'("MAE per unit volume : ",G18.10)') de/omega
  close(71)
  open(50,file='MAE.OUT',form='FORMATTED')
  write(50,'(G18.10)') de
  close(50)
  open(50,file='MAEPUV.OUT',form='FORMATTED')
  write(50,'(G18.10)') de/omega
  close(50)
  write(*,*)
  write(*,'("Info(mae):")')
  write(*,'(" Estimated magnetic anisotropy energy written to MAE.OUT")')
  write(*,'(" MAE per unit volume written to MAEPUV.OUT")')
  write(*,*)
  write(*,'(" Number of fixed spin moment directions used : ",I6)') npmae
  write(*,*)
  write(*,'(" Additional information written to MAE_INFO.OUT")')
end if
! write the MAE to test file
call writetest(28,'magnetic anisotropy energy',tol=1.d-5,rv=de)
! restore original input parameters
avec(:,:)=avec0(:,:)
spinpol=spinpol0
spinorb=spinorb0
cmagz=cmagz0
fsmtype=fsmtype0
bfieldc0(:)=bfieldc00(:)
reducebf=reducebf0
vkloff(:)=vkloff0(:)
trotsht=.false.
end subroutine

