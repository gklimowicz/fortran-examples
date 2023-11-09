
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnit(nmatp,ngp,igpig,vpl,vgpl,vgpc,apwalm,evalfv,evecfv)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: nmatp,ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vgpl(3,ngkmax),vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(out) :: evalfv(nstfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv)
! local variables
integer ns,ist,it,i
integer is,ias,nthd
integer lwork,info
real(8) rmax,t1
real(8) ts1,ts0
! allocatable arrays
real(8), allocatable :: w(:),rwork(:)
complex(8), allocatable :: h(:,:),o(:,:),hv(:,:),ov(:,:)
complex(8), allocatable :: u(:,:),hu(:,:),ou(:,:)
complex(8), allocatable :: hs(:,:),os(:,:),work(:)
! external functions
complex(8), external :: zdotc
ns=2*nstfv
if (iscl.ge.2) then
! read in the eigenvalues/vectors from file
  call getevalfv(filext,0,vpl,evalfv)
  call getevecfv(filext,0,vpl,vgpl,evecfv)
else
! initialise the eigenvectors to canonical basis vectors
  evecfv(1:nmatp,:)=0.d0
  do ist=1,nstfv
    evecfv(ist,ist)=1.d0
  end do
end if
! compute Hamiltonian and overlap matrices
call timesec(ts0)
allocate(h(nmatp,nmatp),o(nmatp,nmatp))
call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP PRIVATE(i,ias,is) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
! Hamiltonian
do i=1,nmatp
  h(1:i,i)=0.d0
end do
do ias=1,natmtot
  is=idxis(ias)
  call hmlaa(tefvr,is,ias,ngp,apwalm(:,:,:,ias),nmatp,h)
  call hmlalo(is,ias,ngp,apwalm(:,:,:,ias),nmatp,h)
  call hmllolo(is,ias,ngp,nmatp,h)
end do
call hmlistl(ngp,igpig,vgpc,nmatp,h)
!$OMP SECTION
! overlap
do i=1,nmatp
  o(1:i,i)=0.d0
end do
do ias=1,natmtot
  is=idxis(ias)
  call olpaa(tefvr,is,ngp,apwalm(:,:,:,ias),nmatp,o)
  call olpalo(is,ias,ngp,apwalm(:,:,:,ias),nmatp,o)
  call olplolo(is,ias,ngp,nmatp,o)
end do
call olpistl(ngp,igpig,nmatp,o)
!$OMP END PARALLEL SECTIONS
call freethd(nthd)
call timesec(ts1)
!$OMP ATOMIC
timemat=timemat+ts1-ts0
call timesec(ts0)
allocate(w(ns),rwork(3*ns))
allocate(hv(nmatp,nstfv),ov(nmatp,nstfv))
allocate(u(nmatp,nstfv),hu(nmatp,nstfv),ou(nmatp,nstfv))
allocate(hs(ns,ns),os(ns,ns))
lwork=2*ns
allocate(work(lwork))
call holdthd(nstfv,nthd)
! start iteration loop
do it=1,maxitefv
  rmax=0.d0
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(t1) &
!$OMP NUM_THREADS(nthd)
  do ist=1,nstfv
! operate with O on the current eigenvector
    call zhemv('U',nmatp,zone,o,nmatp,evecfv(:,ist),1,zzero,ov(:,ist),1)
! normalise the eigenvector
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,ov(:,ist),1))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      call zdscal(nmatp,t1,evecfv(:,ist),1)
      call zdscal(nmatp,t1,ov(:,ist),1)
    end if
! operate with H on the current eigenvector
    call zhemv('U',nmatp,zone,h,nmatp,evecfv(:,ist),1,zzero,hv(:,ist),1)
! estimate the eigenvalue
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,hv(:,ist),1))
    if ((iscl.le.1).and.(it.eq.1)) then
      evalfv(ist)=t1
    else
      evalfv(ist)=(1.d0-befvit)*evalfv(ist)+befvit*t1
    end if
! compute the residual |u> = (H - eO)|v>
    call zcopy(nmatp,hv(:,ist),1,u(:,ist),1)
    t1=-evalfv(ist)
    u(:,ist)=u(:,ist)+t1*ov(:,ist)
! apply the overlap matrix to the residual
    call zhemv('U',nmatp,zone,o,nmatp,u(:,ist),1,zzero,ou(:,ist),1)
! compute the overlap of the residual with itself
    t1=dble(zdotc(nmatp,u(:,ist),1,ou(:,ist),1))
!$OMP ATOMIC
    rmax=max(rmax,t1)
! normalise the residual
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      call zdscal(nmatp,t1,u(:,ist),1)
      call zdscal(nmatp,t1,ou(:,ist),1)
    end if
! apply the Hamiltonian matrix to the residual
    call zhemv('U',nmatp,zone,h,nmatp,u(:,ist),1,zzero,hu(:,ist),1)
  end do
!$OMP END PARALLEL DO
! compute the Hamiltonian and overlap matrices in the subspace formed by the
! eigenvectors and their residuals
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
  do ist=1,nstfv
    call zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,hv(:,ist),1,zzero, &
     hs(1,ist),1)
    call zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,hu(:,ist),1,zzero, &
     hs(1,nstfv+ist),1)
    call zgemv('C',nmatp,nstfv,zone,u,nmatp,hu(:,ist),1,zzero, &
     hs(nstfv+1,nstfv+ist),1)
  end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
  do ist=1,nstfv
    call zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,ov(:,ist),1,zzero, &
     os(1,ist),1)
    call zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,ou(:,ist),1,zzero, &
     os(1,nstfv+ist),1)
    call zgemv('C',nmatp,nstfv,zone,u,nmatp,ou(:,ist),1,zzero, &
     os(nstfv+1,nstfv+ist),1)
  end do
!$OMP END PARALLEL DO
! solve the generalised eigenvalue problem in the subspace
  call zhegv(1,'V','U',ns,hs,ns,os,ns,w,work,lwork,rwork,info)
  if (info.ne.0) exit
! construct the new eigenvectors
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
  do ist=1,nstfv
    call zgemv('N',nmatp,nstfv,zone,evecfv,nmatmax,hs(1,ist),1,zzero, &
     ov(:,ist),1)
    call zgemv('N',nmatp,nstfv,zone,u,nmatp,hs(nstfv+1,ist),1,zone,ov(:,ist),1)
  end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
  do ist=1,nstfv
    call zcopy(nmatp,ov(:,ist),1,evecfv(:,ist),1)
  end do
!$OMP END PARALLEL DO
! check for convergence
  rmax=sqrt(abs(rmax)/dble(nmatp))
  if ((it.ge.minitefv).and.(rmax.lt.epsefvit)) exit
! end iteration loop
end do
call freethd(nthd)
deallocate(w,rwork,h,o,hv,ov)
deallocate(u,hu,ou,hs,os,work)
call timesec(ts1)
!$OMP ATOMIC
timefv=timefv+ts1-ts0
end subroutine

