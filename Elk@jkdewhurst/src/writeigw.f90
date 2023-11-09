
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeigw(fnum)
use modmain
use modgw
implicit none
! arguments
integer fnum
write(fnum,'("+----------------------------------------|")')
write(fnum,'("| Self-consistent density GW calculation |")')
write(fnum,'("+----------------------------------------|")')
write(fnum,*)
write(fnum,'("Temperature (K) : ",G18.10)') tempk
write(fnum,*)
write(fnum,'("Matsubara frequency cut-off : ",G18.10)') wmaxgw
write(fnum,'("Number of Matsubara frequencies : ",I6)') nwgw
write(fnum,*)
write(fnum,'("Maximum |G| for response function : ",G18.10)') gmaxrf
write(fnum,'("Number of response-function G-vectors : ",I8)') ngrf
flush(fnum)
end subroutine

