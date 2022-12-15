
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine maguplot
use modmain
use modulr
use modomp
implicit none
! local variables
integer ifq,idm,nthd
! allocatable arrays
complex(8), allocatable :: magqir_(:,:,:)
! initialise universal variables
call init0
call init1
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(maguplot): spin-unpolarised calculation")')
  write(*,*)
  stop
end if
! initialise the ultra long-range variables
call initulr
! read in the magnetisation from STATE_ULR.OUT
call readstulr
! convert interstitial magnetisation from coarse to fine grid
allocate(magqir_(ngtot,ndmag,nfqrz))
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(idm) &
!$OMP NUM_THREADS(nthd)
do ifq=1,nfqrz
  do idm=1,ndmag
    call zfirctof(magqir(:,idm,ifq),magqir_(:,idm,ifq))
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! write the magnetisation plot to file
select case(task)
case(771)
  open(50,file='MAGU1D.OUT',form='FORMATTED')
  open(51,file='MAGULINES.OUT',form='FORMATTED')
  call plotu1d(50,51,ndmag,magqmt,magqir_)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(maguplot):")')
  write(*,'(" 1D ultra long-range magnetisation plot written to MAGU1D.OUT")')
  write(*,'(" vertex location lines written to MAGULINES.OUT")')
case(772)
  open(50,file='MAGU2D.OUT',form='FORMATTED')
  call plotu2d(.true.,50,ndmag,magqmt,magqir_)
  close(50)
  write(*,*)
  write(*,'("Info(maguplot): 2D ultra long-range magnetisation plot written to &
   &MAGU2D.OUT")')
  if (ndmag.eq.3) then
    write(*,'(" Note that the 3D vector field has been locally projected")')
    write(*,'(" onto the 2D plotting plane axes")')
  end if
case(773)
  open(50,file='MAGU3D.OUT',form='FORMATTED')
  call plotu3d(50,ndmag,magqmt,magqir_)
  close(50)
  write(*,*)
  write(*,'("Info(maguplot): 3D ultra long-range magnetisation plot written to &
   &MAGU3D.OUT")')
end select
deallocate(magqir_)
end subroutine

