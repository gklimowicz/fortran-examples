
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine dysonr(ik,wr,sem,sf)
use modmain
use modgw
use modomp
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: wr(nwplot)
complex(8), intent(in) :: sem(nstsv,nstsv,0:nwfm)
real(8), intent(out) :: sf(nwplot)
! local variables
integer ist,jst,iw
integer nthd
real(8) w,e,sm,t1
complex(8) z1
! allocatable arrays
complex(8), allocatable :: ser(:,:,:),gs(:),g(:,:)
allocate(ser(nstsv,nstsv,nwplot))
ser(:,:,:)=0.d0
call holdthd(nstsv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ist) &
!$OMP NUM_THREADS(nthd)
do jst=1,nstsv
  do ist=1,nstsv
    if (tsediag.and.(ist.ne.jst)) cycle
! perform analytic continuation from the imaginary to the real axis
    call acgwse(ist,jst,sem,wr,ser)
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! solve the Dyson equation for each frequency
call holdthd(nwplot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(gs,g,w,ist,jst) &
!$OMP PRIVATE(e,t1,z1,sm) &
!$OMP NUM_THREADS(nthd)
allocate(gs(nstsv),g(nstsv,nstsv))
!$OMP DO
do iw=1,nwplot
  w=wr(iw)
! compute the diagonal matrix G_s
  do ist=1,nstsv
    e=evalsv(ist,ik)-efermi
    t1=sign(swidth,e)
    gs(ist)=1.d0/cmplx(w-e,t1,8)
  end do
! compute 1 - G_s Sigma
  do ist=1,nstsv
    z1=-gs(ist)
    g(ist,:)=z1*ser(ist,:,iw)
    g(ist,ist)=g(ist,ist)+1.d0
  end do
! invert this matrix
  call zminv(nstsv,g)
! compute G = (1 - G_s Sigma)^(-1) G_s
  do jst=1,nstsv
    z1=gs(jst)
    g(:,jst)=g(:,jst)*z1
  end do
! determine the spectral function
  sm=0.d0
  do ist=1,nstsv
    sm=sm+abs(aimag(g(ist,ist)))
  end do
  sf(iw)=sm*occmax/pi
end do
!$OMP END DO
deallocate(gs,g)
!$OMP END PARALLEL
call freethd(nthd)
deallocate(ser)
end subroutine

