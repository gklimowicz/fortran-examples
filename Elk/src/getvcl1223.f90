
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getvcl1223
! !INTERFACE:
subroutine getvcl1223(ikp,vcl1223)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced set (in,integer)
!   vcl1223 : Coulomb matrix elements (out,complex(nstsv,nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Retrieves Coulomb matrix elements of the type $V(1,2,2,3)$ from the file
!   {\tt VCL1223.OUT}.
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(out) :: vcl1223(nstsv,nstsv,nstsv,nkpt)
! local variables
integer recl
! determine record length
inquire(iolength=recl) vcl1223
!$OMP CRITICAL(u262)
open(262,file='VCL1223.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(262,rec=ikp) vcl1223
close(262)
!$OMP END CRITICAL(u262)
end subroutine
!EOC

