
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initoep
use modmain
implicit none
! local variables
integer is,ic,ist,m
! find maximum core states over all species
ncrmax=0
do is=1,nspecies
  ic=0
  do ist=1,nstsp(is)
    if (spcore(ist,is)) then
      do m=-ksp(ist,is),ksp(ist,is)-1
        ic=ic+1
      end do
    end if
  end do
  ncrmax=max(ncrmax,ic)
end do
! allocate the exchange potential and magnetic field
if (allocated(vxmt)) deallocate(vxmt)
allocate(vxmt(npcmtmax,natmtot))
if (allocated(vxir)) deallocate(vxir)
allocate(vxir(ngtot))
if (spinpol) then
  if (allocated(bxmt)) deallocate(bxmt)
  allocate(bxmt(npcmtmax,natmtot,ndmag))
  if (allocated(bxir)) deallocate(bxir)
  allocate(bxir(ngtot,ndmag))
end if
! zero the exchange potential and magnetic field
vxmt(:,:)=0.d0
vxir(:)=0.d0
if (spinpol) then
  bxmt(:,:,:)=0.d0
  bxir(:,:)=0.d0
end if
end subroutine

