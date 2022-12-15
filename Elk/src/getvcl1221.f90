
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getvcl1221
! !INTERFACE:
subroutine getvcl1221(ikp,vcl1221)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced set (in,integer)
!   vcl1221 : Coulomb matrix elements (out,real(nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Retrieves Coulomb matrix elements of the type $V(1,2,2,1)$ from the file
!   {\tt VCL1221.OUT}.
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(out) :: vcl1221(nstsv,nstsv,nkpt)
! local variables
integer recl
! determine record length
inquire(iolength=recl) vcl1221
!$OMP CRITICAL(u260)
open(260,file='VCL1221.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(260,rec=ikp) vcl1221
close(260)
!$OMP END CRITICAL(u260)
end subroutine
!EOC

