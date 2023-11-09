
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradzvcln(is,gzfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: is
complex(8), intent(out) :: gzfmt(npmtmax,3)
! local variables
integer nr,nri,iro,i0,i1
! automatic arrays
complex(8) zvclmt(npmtmax)
nr=nrmt(is)
nri=nrmti(is)
iro=nri+1
! convert nuclear Coulomb potential to complex spherical harmonics expansion
zvclmt(1:npmt(is))=0.d0
i1=lmmaxi*(nri-1)+1
zvclmt(1:i1:lmmaxi)=vcln(1:nri,is)
i0=i1+lmmaxi
i1=lmmaxo*(nr-iro)+i0
zvclmt(i0:i1:lmmaxo)=vcln(iro:nr,is)
! compute the gradient of the potential
call gradzfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),zvclmt,npmtmax,gzfmt)
end subroutine

