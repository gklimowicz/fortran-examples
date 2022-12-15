
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potuplot
use modmain
use modulr
implicit none
! local variables
integer ifq,ias,is,npc
! allocatable arrays
complex(8), allocatable :: zfmt(:)
! initialise universal variables
call init0
call init1
! initialise the ultra long-range variables
call initulr
! read in the Kohn-Sham potential from STATE_ULR.OUT
call readstulr
! convert potential to spherical harmonics
allocate(zfmt(npcmtmax))
do ifq=1,nfqrz
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
    zfmt(1:npc)=vsqmt(1:npc,ias,ifq)
    call zfsht(nrcmt(is),nrcmti(is),zfmt,vsqmt(:,ias,ifq))
  end do
end do
deallocate(zfmt)
! write the density plot to file
select case(task)
case(741)
  open(50,file='VSU1D.OUT',form='FORMATTED')
  open(51,file='VSULINES.OUT',form='FORMATTED')
  call plotu1d(50,51,1,vsqmt,vsqir)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(potuplot):")')
  write(*,'(" 1D ultra long-range Kohn-Sham potential plot written to &
   &VSU1D.OUT")')
  write(*,'(" vertex location lines written to VSULINES.OUT")')
case(742)
  open(50,file='VSU2D.OUT',form='FORMATTED')
  call plotu2d(.false.,50,1,vsqmt,vsqir)
  open(50)
  write(*,*)
  write(*,'("Info(potuplot):")')
  write(*,'(" 2D ultra long-range Kohn-Sham potential plot written to &
   &VSU2D.OUT")')
case(743)
  open(50,file='VSU3D.OUT',form='FORMATTED')
  call plotu3d(50,1,vsqmt,vsqir)
  close(50)
  write(*,*)
  write(*,'("Info(potuplot):")')
  write(*,'(" 3D ultra long-range Kohn-Sham potential plot written to &
   &VSU3D.OUT")')
end select
end subroutine

