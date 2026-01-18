
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine hmlephb(ik,b)
use modmain
use modphonon
use modbog
use modomp
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(out) :: b(nstsv,nstsv)
! local variables
integer iq,jq,ikq,isym,nthd
integer i1,i2,j1,j2,l1,l2
real(8) vl(3),t0
complex(8) z1,z2
! automatic arrays
complex(4) ephmat(nstsv,nstsv,nbph)
complex(8) x(nbph,nstsv),y(nstsv,nbph)
! prefactor
t0=-2.d0*wqptnr*ephscf(1)**2/dengy
b(:,:)=0.d0
if (.not.anomalous) return
! parallel loop over non-reduced q-points
call holdthd(nqptnr,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ephmat,x,y,jq,vl,isym,ikq) &
!$OMP PRIVATE(i1,i2,j1,j2,l1,l2,z1,z2) &
!$OMP REDUCTION(+:b) &
!$OMP NUM_THREADS(nthd)
do iq=1,nqptnr
! equivalent reduced q-point
  jq=ivqiq(ivq(1,iq),ivq(2,iq),ivq(3,iq))
! k+q-vector in lattice coordinates
  vl(:)=vkl(:,ik)+vql(:,iq)
! index to reduced k+q-vector
  call findkpt(vl,isym,ikq)
! read in the electron-phonon matrix elements from file
  call getephmkq(iq,ik,ephmat)
! perform the contraction
  do i2=1,nstsv
    if (abs(evalsv(i2,ik)-efermi).gt.ecutb) cycle
    do j1=1,nstsv
      do l2=1,nbph
        z1=0.d0
        do j2=1,nstsv
          z1=z1+ephmat(j2,i2,l2)*duv(j1,j2,ikq)
        end do
        x(l2,j1)=z1
      end do
    end do
    do l1=1,nbph
      do j1=1,nstsv
        z1=0.d0
        do l2=1,nbph
          z2=dxx(l2,l1,jq)+dwx(l2,l1,jq)
          z1=z1+z2*x(l2,j1)
        end do
        y(j1,l1)=z1
      end do
    end do
    do i1=1,nstsv
      if (bdiag.and.(i1.ne.i2)) cycle
      if (abs(evalsv(i1,ik)-efermi).gt.ecutb) cycle
      z1=0.d0
      do l1=1,nbph
        do j1=1,nstsv
          z1=z1+conjg(ephmat(j1,i1,l1))*y(j1,l1)
        end do
      end do
      b(i1,i2)=b(i1,i2)+t0*z1
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

