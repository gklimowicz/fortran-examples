
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genjlgpr(ngp,gpc,jlgpr)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: gpc(ngp)
real(8), intent(out) :: jlgpr(njcmax,nspecies,ngp)
! local variables
integer ig,is,n,i,nthd
integer nrc,nrci,irc
real(8) t1,t2
! generate spherical Bessel functions on the coarse radial mesh over all species
call holdthd(ngp,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(t1,t2,is,nrc) &
!$OMP PRIVATE(nrci,n,i,irc) &
!$OMP NUM_THREADS(nthd)
do ig=1,ngp
  t1=gpc(ig)
  do is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    n=lmaxi+1
    i=1
    do irc=1,nrci
      t2=t1*rcmt(irc,is)
      call sbessel(lmaxi,t2,jlgpr(i,is,ig))
      i=i+n
    end do
    n=lmaxo+1
    do irc=nrci+1,nrc
      t2=t1*rcmt(irc,is)
      call sbessel(lmaxo,t2,jlgpr(i,is,ig))
      i=i+n
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

