
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine aceplot
use modmain
use modphonon
use modbog
implicit none
! local variables
integer ik,jk,iq,jq
integer i1,i2,i3,ist,i
real(8) ace,vn,xn
! initialise universal variables
call init0
call init1
call init2
call initeph
!----------------------------------------------------------!
!     plot the fermionic anomalous correlation entropy     !
!----------------------------------------------------------!
open(50,file='FACE3D.OUT',form='FORMATTED',action='WRITE')
write(50,'(3I6," : grid size")') ngridk(:)
do i3=0,ngridk(3)-1
  do i2=0,ngridk(2)-1
    do i1=0,ngridk(1)-1
      ik=ivkiknr(i1,i2,i3)
      jk=ivkik(i1,i2,i3)
      ace=0.d0
      do ist=1,nstsv
        vn=vnorm(ist,jk)
        if ((vn.gt.0.d0).and.(vn.lt.1.d0)) then
          ace=ace+vn*log(vn)+(1.d0-vn)*log(1.d0-vn)
        end if
      end do
      ace=-occmax*ace
      write(50,'(4G18.10)') vkc(:,ik),ace
    end do
  end do
end do
close(50)
!--------------------------------------------------------!
!     plot the bosonic anomalous correlation entropy     !
!--------------------------------------------------------!
open(50,file='BACE3D.OUT',form='FORMATTED',action='WRITE')
write(50,'(3I6," : grid size")') ngridq(:)
do i3=0,ngridq(3)-1
  do i2=0,ngridq(2)-1
    do i1=0,ngridq(1)-1
      iq=ivqiqnr(i1,i2,i3)
      jq=ivqiq(i1,i2,i3)
      ace=0.d0
      do i=1,nbph
        xn=xnorm(i,jq)
        if (xn.gt.0.d0) then
          ace=ace+xn*log(xn)-(1.d0+xn)*log(1.d0+xn)
        end if
      end do
      ace=-ace
      write(50,'(4G18.10)') vqc(:,iq),ace
    end do
  end do
end do
close(50)
write(*,*)
write(*,'("Info(aceplot):")')
write(*,'(" 3D fermionic anomalous correlation entropy written to FACE3D.OUT")')
write(*,'(" 3D bosonic anomalous correlation entropy written to BACE3D.OUT")')
end subroutine

