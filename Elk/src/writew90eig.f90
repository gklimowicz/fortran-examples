
! Copyright (C) 2015 Manh Duc Le, 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin
! and Lars Nordstrom. This file is distributed under the terms of the GNU
! General Public License. See the file COPYING for license details.

subroutine writew90eig
use modmain
use modw90
implicit none
! local variables
integer ik,jk,ist,i
real(8) t1
character(256) fname
fname=trim(seedname)//'.eig'
open(50,file=trim(fname),action='WRITE',form='FORMATTED')
! loop over non-reduced k-points
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
  do i=1,num_bands
    ist=idxw90(i)
    t1=evalsv(ist,jk)-efermi
    write(50,'(2I6,G18.10)') i,ik,t1*ha_ev
  end do
end do
close(50)
write(*,*)
write(*,'("Info(writew90eig): created file ",A)') trim(fname)
end subroutine

