
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine putephmat(iq,ik,ephmat)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq,ik
complex(8), intent(in) :: ephmat(nstsv,nstsv,nbph)
! local variables
integer recl,n
! determine the record length
inquire(iolength=recl) vql(:,iq),vkl(:,ik),nstsv,nbph,ephmat
! record number
n=(iq-1)*nkptnr+ik
!$OMP CRITICAL(u240)
open(240,file='EPHMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
write(240,rec=n) vql(:,iq),vkl(:,ik),nstsv,nbph,ephmat
close(240)
!$OMP END CRITICAL(u240)
end subroutine

