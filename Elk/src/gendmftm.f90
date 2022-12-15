
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma, E. K. U. Gross and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmftm
use modmain
use moddftu
implicit none
! local variables
integer is,ia,ias,i
integer l,k,p,r,t
real(8) t0
! automatic arrays
real(8) wkpr(-lmmaxdm:lmmaxdm)
complex(8) dm(lmmaxdm,2,lmmaxdm,2)
! allocate global array
if (allocated(dmftm)) deallocate(dmftm)
allocate(dmftm(lmmaxdm,2,lmmaxdm,2,natmtot))
! zero the fixed tensor moment density matrices
dmftm(:,:,:,:,:)=0.d0
do i=1,ntmfix
  is=itmfix(1,i)
  if (is.gt.nspecies) then
    write(*,*)
    write(*,'("Error(gendmftm): invalid species number : ",I8)') is
    write(*,*)
    stop
  end if
  ia=itmfix(2,i)
  if (ia.gt.natoms(is)) then
    write(*,*)
    write(*,'("Error(gendmftm): invalid atom number : ",I8)') ia
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  ias=idxas(ia,is)
  l=itmfix(3,i)
  if (l.gt.lmaxdm) then
    write(*,*)
    write(*,'("Error(gendmftm): l > lmaxdm ",2I8)') l,lmaxdm
    write(*,'(" for species ",I4," and atom ",I4)') is,ia
    write(*,*)
    stop
  end if
! generate the 3-index density matrix
  k=itmfix(4,i)
  p=itmfix(5,i)
  r=itmfix(6,i)
  t=itmfix(7,i)
  if (abs(t).gt.lmmaxdm) then
    write(*,*)
    write(*,'("Error(gendmftm): invalid t : ",I8)') t
    write(*,'(" for tensor moment entry ",I3)') i
    write(*,*)
    stop
  end if
! scale factor for conventional normalisation
  t0=sqrt(dble((2*l+1)*2))
  wkpr(:)=0.d0
  wkpr(t)=wkprfix(i)/t0
  call tm3todm(l,k,p,r,lmmaxdm,wkpr,dm)
  dmftm(:,:,:,:,ias)=dmftm(:,:,:,:,ias)+dm(:,:,:,:)
end do
end subroutine

