
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine getephmat(iqp,ikp,ephmat)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iqp,ikp
complex(8), intent(out) :: ephmat(nstsv,nstsv,nbph)
! local variables
integer isym,lspl,iq,ik,iv(3)
integer recl,n,nstsv_,nbph_
real(8) vql_(3),vkl_(3),t1
if (iqp.le.nqpt) then
! q-point is in the reduced set
  iq=iqp
  ik=ikp
else
! q-point is not in the reduced set
  call findqpt(vql(:,iqp),isym,iq)
  lspl=lsplsymc(isym)
  call i3mtv(symlat(:,:,lspl),ivk(:,ikp),iv)
  iv(:)=modulo(iv(:),ngridk(:))
  ik=ivkiknr(iv(1),iv(2),iv(3))
end if
! find the record length
inquire(iolength=recl) vql_,vkl_,nstsv_,nbph_,ephmat
! record number
n=(iq-1)*nkptnr+ik
!$OMP CRITICAL(u240)
open(240,file='EPHMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(240,rec=n) vql_,vkl_,nstsv_,nbph_,ephmat
close(240)
!$OMP END CRITICAL(u240)
t1=abs(vql(1,iq)-vql_(1))+abs(vql(2,iq)-vql_(2))+abs(vql(3,iq)-vql_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getephmat): differing vectors for q-point ",I8)') iq
  write(*,'(" current    : ",3G18.10)') vql(:,iq)
  write(*,'(" EPHMAT.OUT : ",3G18.10)') vql_
  write(*,*)
  stop
end if
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getephmat): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EPHMAT.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getephmat): differing nstsv for q-point, k-point ",2I8)') &
   iq,ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" EPHMAT.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
if (nbph.ne.nbph_) then
  write(*,*)
  write(*,'("Error(getephmat): differing nbph for q-point, k-point ",2I8)') &
   iq,ik
  write(*,'(" current    : ",I8)') nbph
  write(*,'(" EPHMAT.OUT : ",I8)') nbph_
  write(*,*)
  stop
end if
end subroutine

