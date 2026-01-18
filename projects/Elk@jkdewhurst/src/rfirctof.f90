
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfirctof(rfirc,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfirc(ngtc)
real(8), intent(out) :: rfir(ngtot)
! automatic arrays
complex(8) zfftc(ngtc)
! allocatable arrays
complex(8), allocatable :: zfft(:)
allocate(zfft(ngtot))
! Fourier transform function on coarse grid to G-space
zfftc(:)=rfirc(:)
call zfftifc(3,ngdgc,-1,zfftc)
! Fourier transform to fine real-space grid
zfft(:)=0.d0
zfft(igfft(1:ngvc))=zfftc(igfc(1:ngvc))
call zfftifc(3,ngridg,1,zfft)
! output real function
rfir(:)=dble(zfft(:))
deallocate(zfft)
end subroutine

