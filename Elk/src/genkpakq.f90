
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genkpakq
use modmain
use modulr
implicit none
! local variables
integer i1,i2,i3,j1,j2,j3
integer ikpa,ik0,ik
integer n1,iq,ifq,ir
real(8) v(3)
! allocatable arrays
real(8), allocatable :: vkl0(:,:),vkc0(:,:),wkpt0(:)
! find next largest FFT-compatible Q-point grid size (radices 2, 3, 5 and 7)
call nfftifc(4,ngridq(1))
call nfftifc(4,ngridq(2))
call nfftifc(4,ngridq(3))
! total number of Q-points
nqpt=ngridq(1)*ngridq(2)*ngridq(3)
! number of complex FFT elements for real-complex transforms
n1=ngridq(1)/2+1
nfqrz=n1*ngridq(2)*ngridq(3)
! integer grid intervals for the Q-points
intq(1,:)=ngridq(:)/2-ngridq(:)+1
intq(2,:)=ngridq(:)/2
! kappa-point grid should be half the Q-point grid
ngridkpa(:)=(ngridq(:)+1)/2
! number of kappa-points
nkpa=ngridkpa(1)*ngridkpa(2)*ngridkpa(3)
! integer grid intervals for the kappa-points
intkpa(1,:)=ngridkpa(:)/2-ngridkpa(:)+1
intkpa(2,:)=ngridkpa(:)/2
! allocate global Q-point arrays
if (allocated(ivq)) deallocate(ivq)
allocate(ivq(3,nqpt))
if (allocated(ivqiq)) deallocate(ivqiq)
allocate(ivqiq(intq(1,1):intq(2,1),intq(1,2):intq(2,2),intq(1,3):intq(2,3)))
if (allocated(iqfft)) deallocate(iqfft)
allocate(iqfft(nqpt))
if (allocated(ifqrz)) deallocate(ifqrz)
allocate(ifqrz(nqpt))
if (allocated(iqrzf)) deallocate(iqrzf)
allocate(iqrzf(nfqrz))
if (allocated(vql)) deallocate(vql)
allocate(vql(3,nqpt))
if (allocated(vqc)) deallocate(vqc)
allocate(vqc(3,nqpt))
! store the kappa-points as the first nkpa entries in the Q-point arrays
iq=0
do i1=intkpa(1,1),intkpa(2,1)
  do i2=intkpa(1,2),intkpa(2,2)
    do i3=intkpa(1,3),intkpa(2,3)
      iq=iq+1
      ivq(1,iq)=i1
      ivq(2,iq)=i2
      ivq(3,iq)=i3
    end do
  end do
end do
! store the remaining Q-points
do i1=intq(1,1),intq(2,1)
  do i2=intq(1,2),intq(2,2)
    do i3=intq(1,3),intq(2,3)
      if ((i1.lt.intkpa(1,1)).or.(i1.gt.intkpa(2,1)).or. &
          (i2.lt.intkpa(1,2)).or.(i2.gt.intkpa(2,2)).or. &
          (i3.lt.intkpa(1,3)).or.(i3.gt.intkpa(2,3))) then
        iq=iq+1
        ivq(1,iq)=i1
        ivq(2,iq)=i2
        ivq(3,iq)=i3
      end if
    end do
  end do
end do
! ensure the first point is the zero vector
do iq=1,nkpa
  if ((ivq(1,iq).eq.0).and.(ivq(2,iq).eq.0).and.(ivq(3,iq).eq.0)) then
    ivq(:,iq)=ivq(:,1)
    ivq(:,1)=0
    exit
  end if
end do
do iq=1,nqpt
  i1=ivq(1,iq); i2=ivq(2,iq); i3=ivq(3,iq)
! map from (i1,i2,i3) to Q-vector index
  ivqiq(i1,i2,i3)=iq
! Q-vector in Cartesian coordinates
  vqc(:,iq)=dble(i1)*bvecu(:,1) &
           +dble(i2)*bvecu(:,2) &
           +dble(i3)*bvecu(:,3)
! Q-vector in (unit cell) lattice coordinates
  call r3mv(binv,vqc(:,iq),vql(:,iq))
  where(abs(vql(:,iq)).lt.epslat) vql(:,iq)=0.d0
end do
! set up Fourier transform index
do iq=1,nqpt
  i1=ivq(1,iq); i2=ivq(2,iq); i3=ivq(3,iq)
  if (i1.ge.0) then
    j1=i1
  else
    j1=ngridq(1)+i1
  end if
  if (i2.ge.0) then
    j2=i2
  else
    j2=ngridq(2)+i2
  end if
  if (i3.ge.0) then
    j3=i3
  else
    j3=ngridq(3)+i3
  end if
  iqfft(iq)=j3*ngridq(2)*ngridq(1)+j2*ngridq(1)+j1+1
! map from q-point index to real-complex FFT index and vice versa
  if (i1.ge.0) then
    ifq=j3*ngridq(2)*n1+j2*n1+j1+1
    ifqrz(iq)=ifq
    iqrzf(ifq)=iq
  end if
end do
! store the R-vectors in Cartesian coordinates spanning the ultracell
if (allocated(vrcu)) deallocate(vrcu)
allocate(vrcu(3,nqpt))
ir=0
do i3=0,ngridq(3)-1
  v(3)=dble(i3)/dble(ngridq(3))
  do i2=0,ngridq(2)-1
    v(2)=dble(i2)/dble(ngridq(2))
    do i1=0,ngridq(1)-1
      v(1)=dble(i1)/dble(ngridq(1))
      ir=ir+1
      call r3mv(avecu,v,vrcu(:,ir))
    end do
  end do
end do
! store the existing k-point and weight arrays
allocate(vkl0(3,nkpt),vkc0(3,nkpt),wkpt0(nkpt))
vkl0(:,1:nkpt)=vkl(:,1:nkpt)
vkc0(:,1:nkpt)=vkc(:,1:nkpt)
wkpt0(1:nkpt)=wkpt(1:nkpt)
! number of k+kappa-points
nkpt0=nkpt
nkpt=nkpt0*nkpa
! deallocate and reallocate k-point and weight arrays
deallocate(vkl,vkc,wkpt)
allocate(vkl(3,nkpt),vkc(3,nkpt),wkpt(nkpt))
ik=0
do ik0=1,nkpt0
  do ikpa=1,nkpa
    ik=ik+1
    vkl(:,ik)=vkl0(:,ik0)+vql(:,ikpa)
    vkc(:,ik)=vkc0(:,ik0)+vqc(:,ikpa)
    wkpt(ik)=wkpt0(ik0)/dble(nkpa)
  end do
end do
deallocate(vkl0,vkc0,wkpt0)
end subroutine

