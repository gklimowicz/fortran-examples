
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhocoreu
use modmain
use modulr
use modomp
implicit none
! local variables
integer is,ias,ir,i
integer nrc,nrci,irc
integer npc,n,nthd
real(8) t1
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
! generate the core density in spherical coordinates
allocate(rfmt(npcmtmax,natmtot))
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  n=lmmaxi-1
  ir=1
  i=1
  do irc=1,nrci
    t1=rhocr(ir,ias,1)*y00
    rfmt(i:i+n,ias)=t1
    ir=ir+lradstp
    i=i+lmmaxi
  end do
  n=lmmaxo-1
  do irc=nrci+1,nrc
    t1=rhocr(ir,ias,1)*y00
    rfmt(i:i+n,ias)=t1
    ir=ir+lradstp
    i=i+lmmaxo
  end do
end do
! add to the ultra long-range density
call holdthd(nqpt,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ias,is,npc) &
!$OMP NUM_THREADS(nthd)
do ir=1,nqpt
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
    rhormt(1:npc,ias,ir)=rhormt(1:npc,ias,ir)+rfmt(1:npc,ias)
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
deallocate(rfmt)
end subroutine

