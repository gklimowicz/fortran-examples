
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine nesting
use modmain
use modomp
implicit none
! local variables
integer iq,ik,jk,jkq,ivkq(3)
integer ist,i1,i2,i3,nthd
real(8) sm0,sm1,sm2,sm3
real(8) vl(3),vc(3),x,t1
! allocatable arrays
real(8), allocatable :: nq(:)
! external functions
real(8), external :: sdelta
! initialise universal variables
call init0
call init1
call init2
! read Fermi energy from file
call readfermi
! get the eigenvalues from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
end do
allocate(nq(nqpt))
t1=1.d0/swidth
sm0=0.d0
call holdthd(nqpt,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(sm1,sm2,sm3,ik,jk) &
!$OMP PRIVATE(ivkq,jkq,ist,x) &
!$OMP REDUCTION(+:sm0) &
!$OMP NUM_THREADS(nthd)
do iq=1,nqpt
!$OMP CRITICAL(nesting_)
  write(*,'("Info(nesting): ",I6," of ",I6," q-points")') iq,nqpt
!$OMP END CRITICAL(nesting_)
  sm1=0.d0
  do ik=1,nkptnr
    jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
    ivkq(:)=ivk(:,ik)+ivq(:,iq)
    ivkq(:)=mod(ivkq(:),ngridk(:))
    jkq=ivkik(ivkq(1),ivkq(2),ivkq(3))
    sm2=0.d0
    do ist=1,nstsv
      x=(efermi-evalsv(ist,jk))*t1
      sm2=sm2+sdelta(stype,x)*t1
    end do
    sm3=0.d0
    do ist=1,nstsv
      x=(efermi-evalsv(ist,jkq))*t1
      sm3=sm3+sdelta(stype,x)*t1
    end do
    sm1=sm1+sm2*sm3
  end do
  nq(iq)=occmax*omegabz*wkptnr*sm1
  sm0=sm0+omegabz*wqpt(iq)*nq(iq)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
open(50,file='NEST3D.OUT',form='FORMATTED')
write(50,'(3I6," : grid size")') ngridq(:)
do i3=0,ngridq(3)-1
  vl(3)=dble(i3)/dble(ngridq(3))
  do i2=0,ngridq(2)-1
    vl(2)=dble(i2)/dble(ngridq(2))
    do i1=0,ngridq(1)-1
      vl(1)=dble(i1)/dble(ngridq(1))
      vc(:)=bvec(:,1)*vl(1)+bvec(:,2)*vl(2)+bvec(:,3)*vl(3)
      iq=ivqiq(i1,i2,i3)
      write(50,'(4G18.10)') vc(:),nq(iq)
    end do
  end do
end do
close(50)
open(50,file='NESTING.OUT',form='FORMATTED')
write(50,'(G18.10)') sm0
close(50)
write(*,*)
write(*,'("Info(nesting):")')
write(*,'(" Nesting function N(q) written to NEST3D.OUT for plotting")')
write(*,*)
write(*,'(" Total integrated nesting per unit volume written to NESTING.OUT")')
deallocate(nq)
end subroutine

