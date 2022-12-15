
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genstrain
use modmain
implicit none
! local variables
integer i,j,k
real(8) a(3,3),b(3,3),t1
! set first strain equal to isotropic scaling
t1=norm2(avec(:,:))
strain(:,:,1)=avec(:,:)/t1
nstrain=1
do i=1,3
  do j=1,3
! set strain tensor in lattice coordinates to delta_ij
    a(:,:)=0.d0
    a(i,j)=1.d0
! symmetrise strain tensor
    call symmat(a)
! convert to mixed Cartesian-lattice coordinates
    call r3mtm(ainv,a,b)
! orthogonalise strain tensor to previous tensors
    do k=1,nstrain
      t1=sum(b(:,:)*strain(:,:,k))
      b(:,:)=b(:,:)-t1*strain(:,:,k)
    end do
! compute the norm
    t1=norm2(b(:,:))
    if (t1.lt.epslat) cycle
! normalise tensor and store in global array
    nstrain=nstrain+1
    strain(:,:,nstrain)=b(:,:)/t1
  end do
end do
! zero small components
do k=1,nstrain
  do i=1,3
    do j=1,3
      if (abs(strain(i,j,k)).lt.epslat) strain(i,j,k)=0.d0
    end do
  end do
end do
! lattice optimisation case
if ((task.eq.2).or.(task.eq.3)) then
  if (latvopt.eq.2) then
! remove isotropic scaling when latvopt=2
    strain(:,:,1)=strain(:,:,nstrain)
    nstrain=nstrain-1
  else if (latvopt.lt.0) then
! optimise over particular strain when latvopt < 0
    i=abs(latvopt)
    if (i.gt.nstrain) then
      write(*,*)
      write(*,'("Error(genstrain): |latvopt| > nstrain : ",2I8)') i,nstrain
      write(*,*)
      stop
    end if
    strain(:,:,1)=strain(:,:,i)
    nstrain=1
  end if
end if
end subroutine

