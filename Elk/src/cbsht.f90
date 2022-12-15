
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine cbsht(nr,nri,cfmt1,cfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(4), intent(in) :: cfmt1(*)
complex(4), intent(out) :: cfmt2(*)
! local variables
integer i
! transform the inner part of the muffin-tin
call cgemm('N','N',lmmaxi,nri,lmmaxi,cone,cbshti,lmmaxi,cfmt1,lmmaxi,czero, &
 cfmt2,lmmaxi)
! transform the outer part of the muffin-tin
i=lmmaxi*nri+1
call cgemm('N','N',lmmaxo,nr-nri,lmmaxo,cone,cbshto,lmmaxo,cfmt1(i),lmmaxo, &
 czero,cfmt2(i),lmmaxo)
end subroutine

