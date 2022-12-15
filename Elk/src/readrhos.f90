
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readrhos
use modmain
use modtddft
implicit none
integer ios
integer natmtot_,npmtmax_,ngtot_
! allocate static density and charge global arrays
if (allocated(rhosmt)) deallocate(rhosmt)
allocate(rhosmt(npmtmax,natmtot,3))
if (allocated(rhosir)) deallocate(rhosir)
allocate(rhosir(ngtot,3))
if (allocated(chgsmt)) deallocate(chgsmt)
allocate(chgsmt(natmtot,3))
open(100,file='RHOSTAT.OUT',form='UNFORMATTED',action='READ',status='OLD', &
 iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(readrhos): error opening RHOSTAT.OUT")')
  write(*,*)
  stop
end if
read(100) natmtot_
if (natmtot.ne.natmtot_) then
  write(*,*)
  write(*,'("Error(readrhos): differing natmtot")')
  write(*,'(" current     : ",I6)') natmtot
  write(*,'(" RHOSTAT.OUT : ",I6)') natmtot_
  write(*,*)
  stop
end if
read(100) npmtmax_
if (npmtmax.ne.npmtmax_) then
  write(*,*)
  write(*,'("Error(readrhos): differing npmtmax")')
  write(*,'(" current     : ",I6)') npmtmax
  write(*,'(" RHOSTAT.OUT : ",I6)') npmtmax_
  write(*,*)
  stop
end if
read(100) ngtot_
if (ngtot.ne.ngtot_) then
  write(*,*)
  write(*,'("Error(readrhos): differing ngtot")')
  write(*,'(" current     : ",I8)') ngtot
  write(*,'(" RHOSTAT.OUT : ",I8)') ngtot_
  write(*,*)
  stop
end if
read(100) rhosmt,rhosir
read(100) chgsmt
read(100) chgstot
close(100)
end subroutine

