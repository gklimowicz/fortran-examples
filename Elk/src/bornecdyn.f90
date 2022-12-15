
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bornecdyn
use modmain
use modphonon
use modtddft
use modmpi
use modtest
implicit none
! local variables
integer its,iw,i
real(8) vc(3),w1,w2,t0,t1,t2
character(256) fext
! allocatable arrays
real(8), allocatable :: w(:),wt(:),jt(:,:)
real(8), allocatable :: f1(:),f2(:)
complex(8), allocatable :: bec(:,:)
! store original parameters
atposl0(:,:,:)=atposl(:,:,:)
afieldc0(:)=afieldc(:)
! no shifting of atomic basis allowed
tshift=.false.
! initialise universal variables
call init0
atposc0(:,:,:)=atposc(:,:,:)
! generate the time step grid
call gentimes
! allocate local arrays
allocate(w(nwplot),wt(ntimes),jt(3,ntimes))
allocate(f1(ntimes),f2(ntimes),bec(nwplot,3))
! generate energy grid (always non-negative)
w1=max(wplot(1),0.d0)
w2=max(wplot(2),w1)
t1=(w2-w1)/dble(nwplot)
do iw=1,nwplot
  w(iw)=w1+t1*dble(iw-1)
end do
! determine the weights for the spline integration
call wsplint(ntimes,times,wt)
! generate a zero A-field and write to file
npulse=0
nramp=0
if (mp_mpi) call genafieldt
! initial ground-state run should start from atomic densities
trdstate=.false.
! begin new Born effective charge task
10 continue
call bectask(80,fext)
! if nothing more to do then restore original input parameters and return
if (isph.eq.0) then
  filext='.OUT'
  atposl(:,:,:)=atposl0(:,:,:)
  afieldc(:)=afieldc0(:)
  tshift=tshift0
  trddatpos=.false.
  deallocate(w,wt,jt,f1,f2,bec)
  return
end if
if (mp_mpi) then
  write(*,'("Info(bornecdyn): working on ",A)') 'BEC'//trim(fext)
end if
! break the crystal symmetry for the displaced atom
atposl(:,:,:)=atposl0(:,:,:)
atposc(:,:,:)=atposc0(:,:,:)
vc(:)=atposc(:,iaph,isph)
vc(ipph)=vc(ipph)-0.5d0*deltaph
call r3mv(ainv,vc,atposl(:,iaph,isph))
! apply a small static A-field
afieldc(:)=0.d0
afieldc(ipph)=1.d-4
! run the ground-state calculation
call gndstate
! subsequent calculations will read in the previous potential
trdstate=.true.
! write zero atomic forces to file
if (mp_mpi) call becforce
! write displacement to file
datposc(:,:,:,:)=0.d0
datposc(ipph,0,iaph,isph)=deltaph
if (mp_mpi) call writedatposc
trddatpos=.true.
! run the time evolution calculation with Ehrenfest dynamics
task=462
call tddft
task=478
! read in the total current from file
call readjtot(jt)
! filter the high-frequency components from the current
do its=1,ntimes
  t1=exp(-swidth*times(its))
  jt(:,its)=t1*jt(:,its)
end do
! compute the dynamical BEC from the Fourier transformed current
t0=1.d0/(deltaph*cos(tdphi))
do i=1,3
  do iw=1,nwplot
    do its=1,ntimes
      t1=jt(i,its)
      t2=w(iw)*times(its)
      f1(its)=t1*cos(t2)
      f2(its)=t1*sin(t2)
    end do
    t1=dot_product(wt(:),f1(:))
    t2=dot_product(wt(:),f2(:))
    bec(iw,i)=t0*cmplx(t1,t2,8)
  end do
end do
! static and nuclear charge
t1=sum(chgsmt(iasph,1:3))/3.d0+spzn(isph)
! write Born effective charge matrix row to file
if (mp_mpi) then
  do i=1,3
    if (i.eq.ipph) then
      t2=t1
    else
      t2=0.d0
    end if
    do iw=1,nwplot
      write(80,'(2G18.10)') w(iw),dble(bec(iw,i))+t2
    end do
    write(80,*)
    do iw=1,nwplot
      write(80,'(2G18.10)') w(iw),aimag(bec(iw,i))
    end do
    write(80,*)
  end do
  close(80)
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! write test file if required and return
if (test) then
  call writetest(478,'dynamical Born effective charge',nv=nwplot,tol=1.d-2, &
   zva=bec)
  return
end if
goto 10
end subroutine

