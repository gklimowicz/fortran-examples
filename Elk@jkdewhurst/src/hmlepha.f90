
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine hmlepha(ik,a)
use modmain
use modphonon
use modbog
use modomp
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(out) :: a(nstsv,nstsv)
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
a(:,:)=0.d0
if (anomalous) goto 10
! parallel loop over non-reduced q-points
call holdthd(nqptnr,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ephmat,x,y,jq,vl,isym,ikq) &
!$OMP PRIVATE(i1,i2,j1,j2,l1,l2,z1,z2) &
!$OMP REDUCTION(+:a) &
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
    do j2=1,nstsv
      do l2=1,nbph
        z1=0.d0
        do j1=1,nstsv
! swap indices of VV† to get the density matrix at k+q
          z1=z1+ephmat(j1,i2,l2)*dvv(j2,j1,ikq)
        end do
        x(l2,j2)=z1
      end do
    end do
    do l1=1,nbph
      do j2=1,nstsv
        z1=0.d0
        do l2=1,nbph
          z2=dxx(l2,l1,jq)+dwx(l2,l1,jq)
          z1=z1+z2*x(l2,j2)
        end do
        y(j2,l1)=z1
      end do
    end do
    do i1=1,i2
      z1=0.d0
      do l1=1,nbph
        do j2=1,nstsv
          z1=z1+conjg(ephmat(j2,i1,l1))*y(j2,l1)
        end do
      end do
      a(i1,i2)=a(i1,i2)+t0*z1
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
10 continue
! add the second-variational eigenvalues minus the Fermi energy
do i1=1,nstsv
  a(i1,i1)=dble(a(i1,i1))+evalsv(i1,ik)-efermi
end do
end subroutine

