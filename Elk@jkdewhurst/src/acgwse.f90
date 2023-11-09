
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine acgwse(ist,jst,sem,wr,ser)
use modmain
use modgw
implicit none
! arguments
integer ist,jst
complex(8), intent(in) :: sem(nstsv,nstsv,0:nwfm)
real(8), intent(in) :: wr(nwplot)
complex(8), intent(out) :: ser(nstsv,nstsv,nwplot)
! allocatable arrays
complex(8), allocatable :: zm(:),zwr(:),zr(:)
allocate(zm(0:nwfm),zwr(nwplot),zr(nwplot))
zm(:)=sem(ist,jst,:)
zwr(:)=wr(:)
select case(actype)
case(1)
! fit a multipole model
  call acpole(zm,zwr,zr)
case(10)
! stabilised Pade approximant
  call pades(nspade,swidth,nwfm+1,wfm,zm,nwplot,zwr,zr)
case default
  write(*,*)
  write(*,'("Error(acgwse): actype not defined : ",I8)') actype
  write(*,*)
  stop
end select
ser(ist,jst,:)=zr(:)
deallocate(zm,zwr,zr)
end subroutine

