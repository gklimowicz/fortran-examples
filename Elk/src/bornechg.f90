
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bornechg
use modmain
use modphonon
use modmpi
use modtest
implicit none
! local variables
integer ip,l
real(8) vc(3),pvl1(3),pvl2(3)
real(8) becl(3),becc(3),t1
! no shifting of the atomic basis
tshift0=tshift
tshift=.false.
! initialise universal variables
call init0
call init1
! store original parameters
atposl0(:,:,:)=atposl(:,:,:)
atposc0(:,:,:)=atposc(:,:,:)
ngridk0(:)=ngridk(:)
maxscl0=maxscl
! begin new Born effective charge task
10 continue
call bectask(80,filext)
! if nothing more to do then restore original input parameters and return
if (isph.eq.0) then
  filext='.OUT'
  tshift=tshift0
  atposl(:,:,:)=atposl0(:,:,:)
  return
end if
if (mp_mpi) then
  write(*,'("Info(bornechg): working on ",A)') 'BEC'//trim(filext)
end if
! dry run: just generate empty BEC files
if (task.eq.209) goto 10
! apply negative atomic displacement
atposl(:,:,:)=atposl0(:,:,:)
atposc(:,:,:)=atposc0(:,:,:)
vc(:)=atposc(:,iaph,isph)
vc(ipph)=vc(ipph)-0.5d0*deltaph
call r3mv(ainv,vc,atposl(:,iaph,isph))
! initial ground-state run should start from atomic densities
trdstate=.false.
! run the ground-state calculation
call gndstate
! subsequent calculations will read in the previous potential
trdstate=.true.
! compute the first polarisation in lattice coordinates
call polar(pvl1)
! apply positive atomic displacement
atposl(:,:,:)=atposl0(:,:,:)
atposc(:,:,:)=atposc0(:,:,:)
vc(:)=atposc(:,iaph,isph)
vc(ipph)=vc(ipph)+0.5d0*deltaph
call r3mv(ainv,vc,atposl(:,iaph,isph))
! run the ground-state calculation again
call gndstate
! compute the second polarisation
call polar(pvl2)
do l=1,3
! add multiple of 2*pi to bring polarisation vectors into coincidence
  pvl1(l)=modulo(pvl1(l),twopi)
  pvl2(l)=modulo(pvl2(l),twopi)
  t1=pvl1(l)-pvl2(l)
  if (abs(t1-twopi).lt.abs(t1)) then
    pvl1(l)=pvl1(l)-twopi
  else if (abs(t1+twopi).lt.abs(t1)) then
    pvl1(l)=pvl1(l)+twopi
  end if
! calculate the Born effective charge from the difference in polarisations
  t1=wkptnr*occmax*dble(nkspolar*ngridk(l))/(twopi*deltaph)
  becl(l)=t1*(pvl2(l)-pvl1(l))
end do
! convert from lattice to Cartesian coordinates
call r3mv(avec,becl,becc)
! add the core and nuclear charge
becc(ipph)=becc(ipph)+chgcr(isph)+spzn(isph)
! write Born effective charge matrix row to file
if (mp_mpi) then
  do ip=1,3
    write(80,'(G18.10," : ip = ",I4)') becc(ip),ip
  end do
  close(80)
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! write test file if required and return
if (test) then
  call writetest(208,'Born effective charge',nv=3,tol=1.d-3,rva=becc)
  return
end if
goto 10
end subroutine

