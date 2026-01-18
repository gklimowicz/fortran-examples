
! Copyright (C) 2002-2009 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dielectric
! !INTERFACE:
subroutine dielectric
! !USES:
use modmain
use modmpi
use modomp
use modtest
! !DESCRIPTION:
!   Computes the dielectric tensor, optical conductivity and plasma frequency.
!   The formulae are taken from {\it Physica Scripta} {\bf T109}, 170 (2004).
!
! !REVISION HISTORY:
!   Created November 2005 (SS and JKD)
!   Added plasma frequency and intraband contribution (S. Lebegue)
!   Complete rewrite, 2008 (JKD)
!   Fixed problem with plasma frequency, 2009 (Marty Blaber and JKD)
!   Parallelised, 2009 (M. Blaber)
!EOP
!BOC
implicit none
! local variables
logical tssr
integer ik,jk,ist,jst
integer iw,ioc,i,j,nthd
real(8) w1,w2,wplas,x
real(8) ei,ej,eji,t1,t2
complex(8) eta,z1
character(256) fname
! allocatable arrays
real(8), allocatable :: w(:)
complex(8), allocatable :: pmat(:,:,:),sigma(:)
! external functions
real(8), external :: sdelta
! initialise universal variables
call init0
call init1
! read Fermi energy from file
call readfermi
! get the eigenvalues and occupancies from file
call readevalsv
call readoccsv
! allocate local arrays
allocate(w(nwplot),sigma(nwplot))
! generate energy grid (always non-negative)
w1=max(wplot(1),0.d0)
w2=max(wplot(2),w1)
t1=(w2-w1)/dble(nwplot)
do iw=1,nwplot
  w(iw)=w1+t1*dble(iw-1)
end do
! i divided by the complex relaxation time
eta=cmplx(0.d0,swidth,8)
! set flag for scissor operator
if (abs(scissor).gt.1.d-8) then
  tssr=.true.
else
  tssr=.false.
end if
! loop over dielectric tensor components
do ioc=1,noptcomp
  i=optcomp(1,ioc)
  j=optcomp(2,ioc)
  wplas=0.d0
  sigma(:)=0.d0
! parallel loop over non-reduced k-points
  call holdthd(nkptnr,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(pmat,jk,ist,jst) &
!$OMP PRIVATE(ei,ej,eji,z1,t1,x) &
!$OMP REDUCTION(+:wplas,sigma) &
!$OMP NUM_THREADS(nthd)
  allocate(pmat(nstsv,nstsv,3))
!$OMP DO
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(dielectric_)
    write(*,'("Info(dielectric): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(dielectric_)
! equivalent reduced k-point
    jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! read in the momentum matrix elements
    call getpmat(vkl(:,ik),pmat)
! valance states
    do ist=1,nstsv
      ei=evalsv(ist,jk)
! conduction states
      do jst=1,nstsv
        ej=evalsv(jst,jk)
        eji=ej-ei
        z1=pmat(ist,jst,i)*conjg(pmat(ist,jst,j))
! scissor correction
        if (tssr) then
          if ((ei.le.efermi).and.(ej.gt.efermi)) then
            t1=(eji+scissor)/eji
            z1=z1*t1**2
            eji=eji+scissor
          end if
        end if
        if (abs(eji).gt.1.d-8) then
          t1=occsv(ist,jk)*(1.d0-occsv(jst,jk)/occmax)/eji
          sigma(:)=sigma(:)+t1*(z1/(w(:)-eji+eta)+conjg(z1)/(w(:)+eji+eta))
        end if
! add to the plasma frequency
        if (intraband) then
          if (i.eq.j) then
            if (ist.eq.jst) then
              x=(ei-efermi)/swidth
              t1=wkptnr*dble(z1)*sdelta(stype,x)/swidth
              wplas=wplas+t1
            end if
          end if
        end if
      end do
    end do
  end do
!$OMP END DO
  deallocate(pmat)
!$OMP END PARALLEL
  call freethd(nthd)
! multiply response function by prefactor
  z1=zi*wkptnr/omega
  sigma(:)=z1*sigma(:)
! add response function and plasma frequency from each process and redistribute
  if (np_mpi.gt.1) then
    call mpi_allreduce(mpi_in_place,sigma,nwplot,mpi_double_complex,mpi_sum, &
     mpicom,ierror)
    call mpi_allreduce(mpi_in_place,wplas,1,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
  end if
! intraband contribution
  if (intraband) then
    if (i.eq.j) then
      wplas=sqrt(occmax*abs(wplas)*fourpi/omega)
! write the plasma frequency to file
      write(fname,'("PLASMA_",2I1,".OUT")') i,j
      open(50,file=trim(fname),form='FORMATTED')
      write(50,'(G18.10," : plasma frequency")') wplas
      close(50)
! add the intraband contribution to sigma
      t1=wplas**2/fourpi
      do iw=1,nwplot
        sigma(iw)=sigma(iw)+t1/(swidth-zi*w(iw))
      end do
    end if
  end if
! write the optical conductivity to file
  write(fname,'("SIGMA_",2I1,".OUT")') i,j
  open(50,file=trim(fname),form='FORMATTED')
  do iw=1,nwplot
    write(50,'(2G18.10)') w(iw),dble(sigma(iw))
  end do
  write(50,*)
  do iw=1,nwplot
    write(50,'(2G18.10)') w(iw),aimag(sigma(iw))
  end do
  close(50)
! write the dielectric function to file
  write(fname,'("EPSILON_",2I1,".OUT")') i,j
  open(50,file=trim(fname),form='FORMATTED')
  t1=0.d0
  if (i.eq.j) t1=1.d0
  do iw=1,nwplot
    t2=t1-fourpi*aimag(sigma(iw)/(w(iw)+eta))
    write(50,'(2G18.10)') w(iw),t2
  end do
  write(50,*)
  do iw=1,nwplot
    t2=fourpi*dble(sigma(iw)/(w(iw)+eta))
    write(50,'(2G18.10)') w(iw),t2
  end do
  close(50)
! end loop over tensor components
end do
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(dielectric):")')
  write(*,'(" dielectric tensor written to EPSILON_ij.OUT")')
  write(*,'(" optical conductivity written to SIGMA_ij.OUT")')
  if (intraband) then
    write(*,'(" plasma frequency written to PLASMA_ij.OUT")')
  end if
  write(*,'(" for components")')
  do ioc=1,noptcomp
    write(*,'("  i = ",I1,", j = ",I1)') optcomp(1:2,ioc)
  end do
end if
! write sigma to test file if required
call writetest(121,'optical conductivity',nv=nwplot,tol=1.d-2,zva=sigma)
deallocate(w,sigma)
end subroutine
!EOC

