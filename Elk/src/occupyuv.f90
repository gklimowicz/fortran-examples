
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine occupyuv
use modmain
use modbog
use modmpi
implicit none
integer ik,ist
real(8) chg,w,vn
real(8) e0,e1,e,t1
! determine the total charge and fermionic anomalous correlation entropy
chg=0.d0
face=0.d0
do ik=1,nkpt
  w=wkpt(ik)
  do ist=1,nstsv
    vn=vnorm(ist,ik)
    chg=chg+w*vn
    if ((vn.gt.0.d0).and.(vn.lt.1.d0)) then
      face=face+w*(vn*log(vn)+(1.d0-vn)*log(1.d0-vn))
    end if
  end do
end do
chg=occmax*chg
face=-occmax*face
! adjust the Fermi energy
efermi=efermi+tauefm*(chgval-chg)
if (mp_mpi) then
  if (abs(chg-chgval).gt.epschg) then
    write(*,*)
    write(*,'("Warning(occupyuv): incorrect charge : ",2G18.10)') chg,chgval
  end if
end if
! estimate the indirect band gap
e0=-1.d8
e1=1.d8
ikgap(1)=1
ikgap(2)=1
do ist=1,nstsv
  do ik=1,nkpt
    e=evaluv(ist,ik)
    if (vnorm(ist,ik).gt.0.5d0) e=-e
    if (e.le.0.d0) then
      if (e.gt.e0) then
        e0=e
        ikgap(1)=ik
      end if
    else
      if (e.lt.e1) then
        e1=e
        ikgap(2)=ik
      end if
    end if
  end do
end do
bandgap(1)=e1-e0
! estimate the direct band gap
e=1.d8
ikgap(3)=1
do ik=1,nkpt
  e0=-1.d8
  e1=1.d8
  do ist=1,nstsv
    t1=evaluv(ist,ik)
    if (vnorm(ist,ik).gt.0.5d0) t1=-t1
    if (t1.le.0.d0) then
      if (t1.gt.e0) e0=t1
    else
      if (t1.lt.e1) e1=t1
    end if
  end do
  t1=e1-e0
  if (t1.lt.e) then
    e=t1
    ikgap(3)=ik
  end if
end do
bandgap(2)=e
end subroutine

