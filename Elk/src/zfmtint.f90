
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure complex(8) function zfmtint(nr,nri,wr,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr)
complex(8), intent(in) :: zfmt(*)
! local variables
integer iro,i0,i1
i1=lmmaxi*(nri-1)+1
zfmtint=sum(wr(1:nri)*zfmt(1:i1:lmmaxi))
iro=nri+1
i0=i1+lmmaxi
i1=lmmaxo*(nr-iro)+i0
zfmtint=zfmtint+sum(wr(iro:nr)*zfmt(i0:i1:lmmaxo))
zfmtint=fourpi*y00*zfmtint
end function

