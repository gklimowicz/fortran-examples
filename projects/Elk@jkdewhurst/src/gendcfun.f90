
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendcfun
use modmain
use modphonon
implicit none
! local variables
integer ig
real(8) v1,v2,v3,t1,t2
complex(8) z1
v1=atposc(1,iaph,isph)
v2=atposc(2,iaph,isph)
v3=atposc(3,iaph,isph)
do ig=1,ngtot
  t1=vgqc(1,ig)*v1+vgqc(2,ig)*v2+vgqc(3,ig)*v3
  t2=ffacgq(ig,isph)*vgqc(ipph,ig)
  z1=t2*cmplx(sin(t1),cos(t1),8)
  dcfunig(ig)=z1
  dcfunir(igfft(ig))=z1
end do
call zfftifc(3,ngridg,1,dcfunir)
end subroutine

