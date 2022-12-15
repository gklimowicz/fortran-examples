
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine hmlephde(iq,d,e)
use modmain
use modphonon
use modbog
use modomp
implicit none
! arguments
integer, intent(in) :: iq
complex(8), intent(out) :: d(nbph,nbph),e(nbph,nbph)
! local variables
integer ik,jk,ikq,isym,nthd
integer i1,i2,j1,j2,l1,l2
real(8) vl(3),t0
complex(8) z1
! automatic arrays
complex(4) ephmat(nstsv,nstsv,nbph)
complex(8) x(nstsv,nstsv),y(nstsv,nstsv)
! prefactor
t0=-occmax*wkptnr*ephscf(1)**2/dengy
e(:,:)=0.d0
! parallel loop over non-reduced k-points
call holdthd(nkptnr,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ephmat,x,y,jk,vl,isym,ikq) &
!$OMP PRIVATE(l1,l2,j1,j2,i1,i2,z1) &
!$OMP REDUCTION(+:e) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! k+q-vector in lattice coordinates
  vl(:)=vkl(:,ik)+vql(:,iq)
! index to reduced k+q-vector
  call findkpt(vl,isym,ikq)
! read in the electron-phonon matrix elements from file
  call getephmkq(iq,ik,ephmat)
! perform the contraction
  if (anomalous) then
    do l2=1,nbph
      do j1=1,nstsv
        do i2=1,nstsv
          z1=0.d0
          do j2=1,nstsv
            z1=z1+ephmat(j2,i2,l2)*duv(j1,j2,ikq)
          end do
          x(i2,j1)=z1
        end do
      end do
      do i1=1,nstsv
        do j1=1,nstsv
          z1=0.d0
          do i2=1,nstsv
            z1=z1-conjg(duv(i1,i2,jk))*x(i2,j1)
          end do
          y(j1,i1)=z1
        end do
      end do
      do l1=1,nbph
        if (ediag.and.(l1.ne.l2)) cycle
        z1=0.d0
        do i1=1,nstsv
          do j1=1,nstsv
            z1=z1+conjg(ephmat(j1,i1,l1))*y(j1,i1)
          end do
        end do
        e(l1,l2)=e(l1,l2)+t0*z1
      end do
    end do
  else
    do l2=1,nbph
      do j2=1,nstsv
        do i2=1,nstsv
          z1=0.d0
          do j1=1,nstsv
! swap indices of VV† to get the density matrix at k+q
            z1=z1+ephmat(j1,i2,l2)*dvv(j2,j1,ikq)
          end do
          x(i2,j2)=z1
        end do
      end do
      do j2=1,nstsv
        do i1=1,nstsv
          z1=0.d0
          do i2=1,nstsv
! swap indices of VV† to get density matrix at k
            z1=z1+dvv(i2,i1,jk)*x(i2,j2)
          end do
          y(i1,j2)=z1
        end do
      end do
      do l1=1,nbph
        if (ediag.and.(l1.ne.l2)) cycle
        z1=0.d0
        do j2=1,nstsv
          do i1=1,nstsv
            z1=z1+conjg(ephmat(j2,i1,l1))*y(i1,j2)
          end do
        end do
        e(l1,l2)=e(l1,l2)+t0*z1
      end do
    end do
  end if
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! determine the matrix D = D0 or D = D0 + E
if (tephde) then
  d(:,:)=e(:,:)
else
  d(:,:)=0.d0
end if
do l1=1,nbph
  d(l1,l1)=d(l1,l1)+wphq(l1,iq)
end do
end subroutine

