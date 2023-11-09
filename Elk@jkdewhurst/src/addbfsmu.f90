
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine addbfsmu
use modmain
use modulr
implicit none
! local variables
integer idm,is,ias,npc
real(8) t1
! add the global fixed spin moment B-field to the Kohn-Sham field
if ((abs(fsmtype).eq.1).or.(abs(fsmtype).eq.3)) then
  do idm=1,ndmag
    t1=bfsmc(idm)
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      bsqmt(1:npc,ias,idm,1)=bsqmt(1:npc,ias,idm,1)+t1
    end do
    bsqir(1:ngtot,idm,1)=bsqir(1:ngtot,idm,1)+t1*cfunir(1:ngtot)
  end do
end if
! add the muffin-tin fields
if ((abs(fsmtype).eq.2).or.(abs(fsmtype).eq.3)) then
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      t1=bfsmcmt(idm,ias)
      bsqmt(1:npc,ias,idm,1)=bsqmt(1:npc,ias,idm,1)+t1
    end do
  end do
end if
end subroutine

