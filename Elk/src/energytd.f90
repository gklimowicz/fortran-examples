
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine energytd
use modmain
use modtddft
implicit none
! local variables
integer is,ias
real(8) ca,engya
! external functions
real(8), external :: rfinp
! Coulomb potential energy
engyvcl=rfinp(rhomt,rhoir,vclmt,vclir)
! Madelung term
engymad=0.d0
do ias=1,natmtot
  is=idxis(ias)
  engymad=engymad+0.5d0*spzn(is)*(vclmt(1,ias)-vcln(1,is))*y00
end do
! exchange and correlation energy
engyx=rfinp(rhomt,rhoir,exmt,exir)
engyc=rfinp(rhomt,rhoir,ecmt,ecir)
! external vector potential interaction energy
ca=-1.d0/solsc
engya=ca*dot_product(afieldt(:,itimes),jtot(:))
! total energy
engytot=engykn+0.5d0*engyvcl+engymad+engyx+engyc+engya
end subroutine

