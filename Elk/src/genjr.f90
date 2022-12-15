
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genjr
use modmain
use modtddft
implicit none
! local variables
integer is,ias,np,l
real(8) ca,t1
! generate the paramagnetic current
call genjpr
! add the diamagnetic term if a time-dependent A-field is present
if (tafieldt) then
! coupling constant of the external A-field (-1/c)
  ca=-1.d0/solsc
! muffin-tin part
  do l=1,3
    t1=ca*afieldt(l,itimes)
    do ias=1,natmtot
      is=idxis(ias)
      np=npmt(is)
      jrmt(1:np,ias,l)=jrmt(1:np,ias,l)+t1*(rhomt(1:np,ias)-rhosmt(1:np,ias,l))
    end do
  end do
! interstitial part
  do l=1,3
    t1=ca*afieldt(l,itimes)
    jrir(1:ngtot,l)=jrir(1:ngtot,l)+t1*(rhoir(1:ngtot)-rhosir(1:ngtot,l))
  end do
end if
end subroutine

