
! Copyright (C) 2015 Manh Duc Le, 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin
! and Lars Nordstrom. This file is distributed under the terms of the GNU
! General Public License. See the file COPYING for license details.

!BOP
! !ROUTINE: writew90win
! !INTERFACE:
subroutine writew90win
! !USES:
use modmain
use modw90
! !DESCRIPTION:
!   Writes out a template {\tt seedname.win} file with the Wannier90 input
!   parameters. Uses {\tt wannier} and {\tt wannierExtra} blocks.
!
! !REVISION HISTORY:
!   Created January 2015 (Manh Duc Le)
!   Modified, August 2018 (Arsenii Gerasimov)
!   Modified, February 2019 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,is,ia,i
character(256) fname
fname=trim(seedname)//'.win'
open(50,file=trim(fname),action='WRITE',form='FORMATTED')
! write the number of Wannier functions and bands
write(50,'("length_unit = bohr")')
write(50,'("num_wann =  ",I8)') num_wann
write(50,'("num_bands = ",I8)') num_bands
write(50,'("num_iter = ",I8)') num_iter
write(50,*)
write(50,'("trial_step = ",G18.10)') trial_step
write(50,*)
write(50,'("use_bloch_phases = true")')
! write spinors if required
if (spinpol) then
  write(50,*)
  write(50,'("spinors = true")')
  write(50,'("spn_formatted = true")')
end if
! write lattice vectors
write(50,*)
write(50,'("begin unit_cell_cart")')
write(50,'("bohr")')
write(50,'(3G18.10)') avec(:,1)
write(50,'(3G18.10)') avec(:,2)
write(50,'(3G18.10)') avec(:,3)
write(50,'("end unit_cell_cart")')
! writes atomic positions
write(50,*)
write(50,'("begin atoms_frac")')
do is=1,nspecies
  do ia=1,natoms(is)
    write(50,'(A5,3G18.10)') trim(spsymb(is)),atposl(:,ia,is)
  end do
end do
write(50,'("end atoms_frac")')
! write the list of k-points
write(50,*)
write(50,'("mp_grid = ",3I6)') ngridk
write(50,*)
write(50,'("begin kpoints")')
do ik=1,nkptnr
  write(50,'(3G18.10)') vkl(:,ik)
end do
write(50,'("end kpoints")')
write(50,*)
! write the extra lines
do i=1,nxlwin
  write(50,'(A)') trim(xlwin(i))
end do
close(50)
write(*,*)
write(*,'("Info(writew90win): created file ",A)') trim(fname)
end subroutine
!EOC

