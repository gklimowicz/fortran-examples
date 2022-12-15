
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writedvs(fext)
use modmain
use modphonon
implicit none
! arguments
character(*), intent(in) :: fext
! local variables
integer is,ias
! allocatable arrays
complex(8), allocatable :: zfmt(:,:,:)
allocate(zfmt(lmmaxo,nrmtmax,natmtot))
open(150,file='DVS'//trim(fext),form='UNFORMATTED',action='WRITE')
write(150) version
write(150) nspecies
write(150) lmmaxo
do is=1,nspecies
  write(150) natoms(is)
  write(150) nrmt(is)
end do
write(150) ngridg
do ias=1,natmtot
  is=idxis(ias)
  call zfmtpack(.false.,nrmt(is),nrmti(is),dvsmt(:,ias),zfmt(:,:,ias))
end do
write(150) zfmt,dvsir
close(150)
deallocate(zfmt)
end subroutine

