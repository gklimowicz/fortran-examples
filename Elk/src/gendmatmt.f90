
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmatmt
use modmain
use moddftu
implicit none
! generate the density matrix in each muffin-tin
call gendmat(.false.,.false.,0,lmaxdm,lmmaxdm,dmatmt)
! initialise with symmetry-breaking tensor moments
if (ftmtype.lt.0) then
  dmftm=dmftm*reducebf
  dmatmt=dmatmt+dmftm
endif
! symmetrise the density matrix
call symdmat(lmaxdm,lmmaxdm,dmatmt)
end subroutine

