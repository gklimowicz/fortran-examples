
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine addbfsm
use modmain
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
      bsmt(1:npc,ias,idm)=bsmt(1:npc,ias,idm)+t1
    end do
    bsir(1:ngtot,idm)=bsir(1:ngtot,idm)+t1*cfunir(1:ngtot)
  end do
end if
! add the muffin-tin fields
if ((abs(fsmtype).eq.2).or.(abs(fsmtype).eq.3)) then
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      t1=bfsmcmt(idm,ias)
      bsmt(1:npc,ias,idm)=bsmt(1:npc,ias,idm)+t1
    end do
  end do
end if
end subroutine

