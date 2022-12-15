
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writestulr
use modmain
use modulr
implicit none
! local variables
integer ifq,idm
open(100,file='STATE_ULR.OUT',form='UNFORMATTED',action='WRITE')
write(100) version
write(100) natmtot
write(100) npcmtmax
write(100) ngtc
write(100) ngtot
write(100) ndmag
write(100) fsmtype
write(100) nqpt
write(100) nfqrz
write(100) ivq
write(100) iqrzf
! write the ultra long-range density in Q-space
do ifq=1,nfqrz
  write(100) rhoqmt(:,:,ifq)
  write(100) rhoqir(:,ifq)
end do
! write the Kohn-Sham effective potential in Q-space
do ifq=1,nfqrz
  write(100) vsqmt(:,:,ifq)
  write(100) vsqir(:,ifq)
end do
! write the external Coulomb potential in Q-space
do ifq=1,nfqrz
  write(100) vclq(ifq)
end do
if (spinpol) then
! write the magnetisation in Q-space
  do ifq=1,nfqrz
    do idm=1,ndmag
      write(100) magqmt(:,:,idm,ifq)
      write(100) magqir(:,idm,ifq)
    end do
  end do
! write the Kohn-Sham effective magnetic field in Q-space
  do ifq=1,nfqrz
    do idm=1,ndmag
      write(100) bsqmt(:,:,idm,ifq)
      write(100) bsqir(:,idm,ifq)
    end do
  end do
! write the external magnetic fields in Q-space
  do ifq=1,nfqrz
    do idm=1,ndmag
      write(100) bfcq(idm,ifq)
      write(100) bfcmtq(:,idm,ifq)
    end do
  end do
! write fixed spin moment magnetic fields
  if (fsmtype.ne.0) then
    write(100) bfsmc
    write(100) bfsmcmt
  end if
end if
close(100)
end subroutine

