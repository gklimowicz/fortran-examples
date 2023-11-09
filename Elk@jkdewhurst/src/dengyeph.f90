
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine dengyeph
use modmain
use modphonon
use modbog
implicit none
! local variables
integer ik,iq,i
real(8) w
! change in electron energy per unit cell
dengye=0.d0
do ik=1,nkpt
  w=wkpt(ik)
  do i=1,nstsv
    dengye=dengye+w*abs((evalsv(i,ik)-efermi)*dble(dvv(i,i,ik)))
  end do
end do
dengye=abs(occmax*dengye)
! change in phonon energy per unit cell
dengyph=0.d0
do iq=1,nqpt
  w=wqpt(iq)
  do i=1,nbph
    dengyph=dengyph+w*abs(wphq(i,iq)*dble(dxx(i,i,iq)))
  end do
end do
! sum of both
dengy=dengye+dengyph
end subroutine

