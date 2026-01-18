
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine ephdos
use modmain
use modphonon
use modbog
implicit none
! local variables
integer nsk(3),ik,ist,iw
real(8) v(3),dw,vn,t1
! allocatable arrays
integer, allocatable :: idx(:)
real(8), allocatable :: w(:),f(:,:),g(:)
! initialise universal variables
call init0
call init1
call init2
call initeph
allocate(idx(nstsv),w(nstsv))
do ik=1,nkpt
! get the eigenvalues from file
  call getevaluv(ik,evaluv(:,ik))
! negate eigenvalues corresponding to V-norm > 1/2
  do ist=1,nstsv
    if (vnorm(ist,ik).gt.0.5d0) evaluv(ist,ik)=-evaluv(ist,ik)
  end do
! arrange in ascending order
  w(:)=evaluv(:,ik)
  call sortidx(nstsv,w,idx)
  evaluv(:,ik)=w(idx(:))
! put the V-norm into the same order as the eigenvalues
  w(:)=vnorm(:,ik)
  vnorm(:,ik)=w(idx(:))
end do
deallocate(idx,w)
! generate the partial and total DOS and write to file
allocate(w(nwplot),f(nstsv,nkptnr),g(nwplot))
! generate frequency grid
dw=(wplot(2)-wplot(1))/dble(nwplot)
do iw=1,nwplot
  w(iw)=dw*dble(iw-1)+wplot(1)
end do
! number of subdivisions used for interpolation in the Brillouin zone
nsk(:)=max(ngrkf/ngridk(:),1)
! set the weight array
f(:,:)=occmax
! integrate over the Brillouin zone
call brzint(nswplot,ngridk,nsk,ivkik,nwplot,wplot,nstsv,nstsv,evaluv,f,g)
! output the total electronic DOS to file
open(50,file='TDOS_EPH.OUT.OUT',form='FORMATTED',action='WRITE')
do iw=1,nwplot
  write(50,'(2G18.10)') w(iw),g(iw)
end do
close(50)
! output the FACE vs energy histogram to file
open(50,file='FACEEH.OUT',form='FORMATTED',action='WRITE')
do ik=1,nkpt
! map k-vector to first Brillouin zone
  v(:)=vkc(:,ik)
  call vecfbz(epslat,bvec,v)
  do ist=1,nstsv
    vn=vnorm(ist,ik)
    if ((vn.gt.0.d0).and.(vn.lt.1.d0)) then
      t1=-(vn*log(vn)+(1.d0-vn)*log(1.d0-vn))
    else
      t1=0.d0
    end if
    if (t1.lt.1.d-4) cycle
    write(50,'(5G18.10)') evaluv(ist,ik),t1,v
  end do
end do
close(50)
write(*,*)
write(*,'("Info(ephdos):")')
write(*,'(" Total electronic density of states for the electron-phonon")')
write(*,'(" system written to TDOS_EPH.OUT.OUT")')
write(*,*)
write(*,'(" Fermionic anomalous correlation entropy vs energy histogram")')
write(*,'(" written to FACEEH.OUT")')
write(*,*)
write(*,'(" Fermi energy is at zero in plots")')
write(*,*)
write(*,'(" DOS units are states/Hartree/unit cell")')
deallocate(w,f,g)
end subroutine

