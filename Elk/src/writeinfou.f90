
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeinfou(fnum)
use modmain
use modulr
implicit none
! arguments
integer fnum
! local variables
real(8) t1
write(fnum,'("+------------------------------+")')
write(fnum,'("| Ultra long-range calculation |")')
write(fnum,'("+------------------------------+")')
write(fnum,*)
write(fnum,'("Ultracell lattice vectors :")')
write(fnum,'(3G18.10)') avecu(1,1),avecu(2,1),avecu(3,1)
write(fnum,'(3G18.10)') avecu(1,2),avecu(2,2),avecu(3,2)
write(fnum,'(3G18.10)') avecu(1,3),avecu(2,3),avecu(3,3)
write(fnum,*)
write(fnum,'("Ultracell reciprocal lattice vectors :")')
write(fnum,'(3G18.10)') bvecu(1,1),bvecu(2,1),bvecu(3,1)
write(fnum,'(3G18.10)') bvecu(1,2),bvecu(2,2),bvecu(3,2)
write(fnum,'(3G18.10)') bvecu(1,3),bvecu(2,3),bvecu(3,3)
write(fnum,*)
write(fnum,'("Ultracell volume                : ",G18.10)') omegau
write(fnum,'("Ultracell Brillouin zone volume : ",G18.10)') omegabzu
write(fnum,*)
t1=omegau/omega
write(fnum,'("Ratio of ultracell to unit cell volume : ",G18.10)') t1
t1=t1*dble(natmtot)
write(fnum,'("Number of atoms in ultracell : ",I16)') nint(t1,8)
write(fnum,*)
write(fnum,'("kappa-point grid : ",3I6)') ngridkpa
write(fnum,'("Q-point grid     : ",3I6)') ngridq
write(fnum,*)
write(fnum,'("Small Q-vector cut-off : ",G18.10)') q0cut
if (fsmtype.ne.0) then
  write(fnum,*)
  write(fnum,'("Fixed spin moment (FSM) calculation, type : ",I4)') fsmtype
  if (fsmtype.lt.0) then
    write(fnum,'("  only moment direction is fixed")')
  end if
end if
write(fnum,*)
write(fnum,'("Hamiltonian matrix size : ",I8)') nstulr
flush(fnum)
end subroutine

