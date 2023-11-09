
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhosplot
use modmain
use modtddft
implicit none
! local variables
integer is,ias,ispn
integer nr,nri,iro,ir,i
! determine the static density and charge
call rhostatic
! remove the core charge density
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
  do ispn=1,nspncr
    i=1
    do ir=1,nri
      rhosmt(i,ias,:)=rhosmt(i,ias,:)-rhocr(ir,ias,ispn)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      rhosmt(i,ias,:)=rhosmt(i,ias,:)-rhocr(ir,ias,ispn)
      i=i+lmmaxo
    end do
  end do
end do
! produce 1D plot of the static density
open(50,file='RHOS1D.OUT',form='FORMATTED',action='WRITE')
open(51,file='RHOSLINES.OUT',form='FORMATTED',action='WRITE')
call plot1d(50,51,3,rhosmt,rhosir)
close(50)
close(51)
write(*,*)
write(*,'("Info(rhosplot):")')
write(*,'(" 1D static density plot written to RHOS1D.OUT")')
write(*,'(" vertex location lines written to RHOSLINES.OUT")')
write(*,*)
write(*,'(" The core density is not included in the plot")')
end subroutine

