
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine writemom(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,ias
write(fnum,*)
write(fnum,'("Moments :")')
write(fnum,'(" interstitial",T30,": ",3G18.10)') momir(1:ndmag)
write(fnum,'(" muffin-tins")')
do is=1,nspecies
  write(fnum,'("  species : ",I4," (",A,")")') is,trim(spsymb(is))
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(fnum,'("   atom ",I4,T30,": ",3G18.10)') ia,mommt(1:ndmag,ias)
  end do
end do
write(fnum,'(" total in muffin-tins",T30,": ",3G18.10)') mommttot(1:ndmag)
write(fnum,'(" total moment",T30,": ",3G18.10)') momtot(1:ndmag)
flush(fnum)
end subroutine

