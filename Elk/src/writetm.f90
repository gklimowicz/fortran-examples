
! Copyright (C) 2009 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetm
use modmain
use moddftu
implicit none
if (dftu.eq.0) then
  write(*,*)
  write(*,'("Error(writetmdu): dftu = 0")')
  write(*,*)
  stop
end if
! initialize universal variables
call init0
! read density matrix from file DMATMT.OUT
call readdmatmt
! generate the DFT+U muffin-tin potential matrices
call genvmatmt
! write tensor moments to TENSMOM.OUT file
call writetm3
write(*,*)
write(*,'("Info(writetm): Tensor moment decomposition of density matrix")')
write(*,'(" in the spherical basis written to TENSMOM.OUT")')
end subroutine

