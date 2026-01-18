
! Copyright (C) 2015 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetddos(fext)
use modmain
use modtddft
implicit none
! arguments
character(*), intent(in) :: fext
! local variables
integer ik,ist,jst
real(8) sm,t1
complex(8) z1
! allocatable arrays
real(8), allocatable :: occsvp(:,:)
complex(8), allocatable :: evecsv(:,:),evecsvt(:,:)
! external functions
complex(8), external :: zdotc
allocate(occsvp(nstsv,nkpt))
do ik=1,nkpt
  allocate(evecsv(nstsv,nstsv),evecsvt(nstsv,nstsv))
! read in ground-state eigenvectors
  call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! read in the time evolving eigenvectors
  call getevecsv('_TD.OUT',ik,vkl(:,ik),evecsvt)
! determine the time-dependent projected occupancies
  do ist=1,nstsv
    sm=0.d0
    do jst=1,nstsv
      t1=occsv(jst,ik)
      if (abs(t1).lt.epsocc) cycle
      z1=zdotc(nstsv,evecsv(:,ist),1,evecsvt(:,jst),1)
      sm=sm+t1*(dble(z1)**2+aimag(z1)**2)
    end do
    occsvp(ist,ik)=sm
  end do
! write projected occupancies to file
  call putoccsv('P'//trim(fext),ik,occsvp(:,ik))
  deallocate(evecsv,evecsvt)
end do
! compute the effective electronic temperature
call tdtemp(occsvp)
! write the DOS to file
call dos(fext,.true.,occsvp)
deallocate(occsvp)
end subroutine

