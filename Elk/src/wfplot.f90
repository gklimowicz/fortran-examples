
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wfplot
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik,ist
real(8) x,t1
! external functions
real(8), external :: sdelta
! initialise universal variables
call init0
call init1
! read the density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! set the occupancies
if ((task.eq.61).or.(task.eq.62).or.(task.eq.63)) then
! plot of wavefunction modulus squared
  ik=kstlist(1,1)
  ist=kstlist(2,1)
  if ((ik.lt.1).or.(ik.gt.nkpt)) then
    write(*,*)
    write(*,'("Error(wfplot): k-point out of range : ",I8)') ik
    write(*,*)
    stop
  end if
  if ((ist.lt.1).or.(ist.gt.nstsv)) then
    write(*,*)
    write(*,'("Error(wfplot): state out of range : ",I8)') ist
    write(*,*)
    stop
  end if
! select a particular wavefunction using its occupancy
  occsv(:,:)=0.d0
  occsv(ist,ik)=1.d0/wkpt(ik)
! no symmetrisation required
  nsymcrys=1
else
! plotting an STM image by setting occupancies to be a delta function at the
! Fermi energy
  t1=1.d0/swidth
  do ik=1,nkpt
! get the eigenvalues from file
    call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
    do ist=1,nstsv
      x=(efermi-evalsv(ist,ik))*t1
      occsv(ist,ik)=occmax*wkpt(ik)*sdelta(stype,x)*t1
    end do
  end do
end if
! compute the valence charge density with the new occupancies
call rhomagv
! write the wavefunction modulus squared plot to file
if (mp_mpi) then
  select case(task)
  case(61)
    open(50,file='WF1D.OUT',form='FORMATTED')
    open(51,file='WFLINES.OUT',form='FORMATTED')
    call plot1d(50,51,1,rhomt,rhoir)
    close(50)
    close(51)
    write(*,*)
    write(*,'("Info(wfplot):")')
    write(*,'(" 1D wavefunction modulus squared written to WF1D.OUT")')
    write(*,'(" vertex location lines written to WFLINES.OUT")')
  case(62)
    open(50,file='WF2D.OUT',form='FORMATTED')
    call plot2d(.false.,50,1,rhomt,rhoir)
    close(50)
    write(*,*)
    write(*,'("Info(wfplot):")')
    write(*,'(" 2D wavefunction modulus squared written to WF2D.OUT")')
  case(162)
    open(50,file='STM2D.OUT',form='FORMATTED')
    call plot2d(.false.,50,1,rhomt,rhoir)
    close(50)
    write(*,*)
    write(*,'("Info(wfplot):")')
    write(*,'(" 2D STM image written to STM2D.OUT")')
  case(63)
    open(50,file='WF3D.OUT',form='FORMATTED')
    call plot3d(50,1,rhomt,rhoir)
    close(50)
    write(*,*)
    write(*,'("Info(wfplot):")')
    write(*,'(" 3D wavefunction modulus squared written to WF3D.OUT")')
  end select
  if (task.ne.162) then
    write(*,'(" for k-point ",I8," and state ",I6)') kstlist(1,1),kstlist(2,1)
  end if
end if
end subroutine

